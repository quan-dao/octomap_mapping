#ifndef OCTOMAP_EVID_OCTREE_H
#define OCTOMAP_EVID_OCTREE_H

#include <octomap/OcTreeNode.h>
#include <octomap/OccupancyOcTreeBase.h>
#include <octomap/octomap_utils.h>

#include <iostream>
#include <vector>

#include <ros/ros.h>


const float EPSILON = 1.0e-5;

namespace octomap {

	// ---------------------- Evidential Mass ---------------------------
	class EvidMass {
	public:
		EvidMass() {
			mass = new float[3];
			mass[0] = 0.0f;  // free
			mass[1] = 0.0f;  // occupied
			mass[2] = 1.0f;  // ignorance 
		}

		EvidMass(float f, float o, float i) {
			mass = new float[3];
			mass[0] = f;  // free
			mass[1] = o;  // occupied
			mass[2] = i;  // ignorance
		}

		~EvidMass()
		{
			delete mass;
		}

		/// accessing functions
		inline float c() const {
			float c_ = 1.0f - (mass[0] + mass[1] + mass[2]);
			if (c_<-EPSILON) // mass leaking because of clamping???
				return 0.0f;
			else
				return c_;
		}
		inline float f() const {return mass[0];}
		inline float o() const {return mass[1];}
		inline float i() const {return mass[2];}
		
		// Check if evidential mass is valid
		inline bool isValid() const {
			bool re = (mass[0]>=-EPSILON) && (mass[1]>=-EPSILON) && (mass[2]>=-EPSILON);
			if (re)
				return true;
			else {
				std::cout<<"[ERROR]Invalid mass ( "<<mass[0]<<", "<<mass[1]<<", "<<mass[2]<<" )\n";
				return false;
			}

		}

		inline void setValue(float f, float o, float i) {
			mass[0] = f;
			mass[1] = o;
			mass[2] = i;
		}

	protected:
		float* mass;  // just store free, occupied & ignorance
	}; 


	// ------------------------ Node definition ------------------------
	// forward declaration for "friend"
	class EvidOcTree;

	class EvidOcTreeNode : public OcTreeNode {
	public:
		friend class EvidOcTree;

		EvidOcTreeNode() : updatedTime(0, 0) {
			massPtr = new EvidMass;
		}

		~EvidOcTreeNode() {
			delete massPtr;
		}

		void copyData(const EvidOcTreeNode& from) {
			// copy data payload
			OcTreeNode::copyData(from);
			// copy mass
			copyMass(from.massPtr);
			// copy updated time
			updatedTime = from.updatedTime;
		}
		
		/// Mass handler
		inline void copyMass(const EvidMass* srcMassPtr) {
			massPtr->setValue(srcMassPtr->f(), srcMassPtr->o(), srcMassPtr->i());
		}

		inline void decayMass(float alp) {
			massPtr->setValue(massPtr->f() * alp,
												massPtr->o() * alp,
												massPtr->i() * alp + 1.0f - alp);
		}

		inline float evidMassToLogOdds() const {
			// use Pignistic transformation to convert evidential mass to occupancy prob
			return logodds((double) (massPtr->o() + 0.5f*massPtr->i()));
		}

		/// Time handler
		inline ros::Time getUpdatedTime() const {return updatedTime;}
		inline void setUpdatedTime(ros::Time& _time) {updatedTime = _time;}

		/**
		 * Update a node's occupancy based on its children's evidential mass
		 * by first computing node's evidential mass as its children's average mass
		 * then convert evidential mass to log odds
		 */
		void updateOccupancyChildren();

		EvidMass* massPtr;
	protected:
		ros::Time updatedTime;
	};
	 

	// ------------------------ Tree definition ------------------------
	class EvidOcTree : public OccupancyOcTreeBase<EvidOcTreeNode> {
	public:
		/// Default constructor: set resolution of leaves
		EvidOcTree(double resolution);
		EvidOcTree(double resolution, float lambda_free_, float lambda_occupied_, float conf_thres_, float tau_);
		/// Virtual constructor creates a new object of the same type
		EvidOcTree* create() const {return new EvidOcTree(resolution);}

		std::string getTreeType() const {return "EvidOcTree";}

		/** Redefinition of updateNode() for it to call functions implementing 
		 * Evidential Fusion instead of Occupancy Fusion
		 */ 
		EvidOcTreeNode* updateNode(const OcTreeKey& key, bool occupied, bool lazy_eval=false);

		/// Alternative for updateNode with log_odds_update
		EvidOcTreeNode* updateNode(const OcTreeKey& key, const EvidMass& basicBeliefAssign, bool lazy_eval=false);

		/// Alternative for updateNodeRecurs with log_odds_update
		EvidOcTreeNode* updateNodeRecurs(EvidOcTreeNode* node, bool node_just_created, const OcTreeKey& key,
                                    unsigned int depth, const EvidMass& basicBeliefAssign, bool lazy_eval=false);

		/**
		 * Counter part of updateNodeLogOdds in Evidential Grid
		 * This function perform 2-step evidential fusion (conjunctive & dempster normalization)
		 * Node's log odds is preserved for pruning & occupancy queries but no longer directly
		 * updated. Instead, node's occupancy probability is retrieved by the Pignistic Transformation
		 */
		void upadteNodeEvidMass(EvidOcTreeNode* node, const EvidMass& basicBeliefAssign, const OcTreeKey& key);
		
		inline void setIncomingTime(const ros::Time& _time) {incomingTime = _time;}
		
		/// conflict_cells_center handle
		inline std::vector<point3d>& getConflictCells() {return conflict_cells_center;} 
		inline void clearConfictCells() {conflict_cells_center.clear();}

	protected:
		// Evidential Fusion constants
		const float tau;  // time constant
		const float lambda_occupied;
		const float lambda_free;
		const float conflict_thres;  // 0.35f
		std::vector<point3d> conflict_cells_center;
		EvidMass basic_belief_occupied, basic_belief_free;

		// timestamp of incoming pointcloud. This is updated in callback function "insertCloudCallback"
		ros::Time incomingTime;  

		/**
     * Static member object which ensures that this OcTree's prototype
     * ends up in the classIDMapping only once. You need this as a 
     * static member in any derived octree class in order to read .ot
     * files through the AbstractOcTree factory. You should also call
     * ensureLinking() once from the constructor.
     */
    class StaticMemberInitializer{
		public:
			StaticMemberInitializer() {
				EvidOcTree* tree = new EvidOcTree(0.05); // 0.1
				tree->clearKeyRays();
				AbstractOcTree::registerTreeType(tree);
			}

			/**
		 * Dummy function to ensure that MSVC does not drop the
		 * StaticMemberInitializer, causing this tree failing to register.
		 * Needs to be called from the constructor of this octree.
		 */
			void ensureLinking() {};
    };
    /// static member to ensure static initialization (only once)
    static StaticMemberInitializer evidOcTreeMemberInit;


	}; 

} //end namespace

#endif
