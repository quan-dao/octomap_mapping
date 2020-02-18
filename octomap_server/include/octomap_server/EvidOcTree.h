#ifndef OCTOMAP_EVID_OCTREE_H
#define OCTOMAP_EVID_OCTREE_H

#include <octomap/OcTreeNode.h>
#include <octomap/OccupancyOcTreeBase.h>
#include <octomap/octomap_utils.h>

#include <iostream>
#include <vector>

#include <ros/ros.h>

namespace octomap {

	// ---------------------- Evidential Mass ---------------------------
	class EvidMass {
	public:
		EvidMass() {
			mass.push_back(0.0f);  // conflict
			mass.push_back(0.0f);  // free
			mass.push_back(0.0f);  // occupied
			mass.push_back(1.0f);  // ignorance 
		}

		EvidMass(float c, float f, float o, float i) {
			mass.push_back(c);
			mass.push_back(f);
			mass.push_back(o);
			mass.push_back(i);
		}

		/// accessing functions
		inline float c() const {return mass[0];}
		inline float f() const {return mass[1];}
		inline float o() const {return mass[2];}
		inline float i() const {return mass[3];}
		
		// Check if evidential mass is valid
		inline bool isValid() const {
			float tot_mass(0.0);
			bool re = true; 
			for(float _m : mass) {
				tot_mass += _m;
				re = re && (_m >= 0.0f);
			}
			return re && (tot_mass >= 0.95f) && (tot_mass <= 1.05f);
		}

		inline void setValue(float c, float f, float o, float i) {
			mass[0] = c;
			mass[1] = f;
			mass[2] = o;
			mass[3] = i;
		}

	protected:
		std::vector<float> mass;
	}; 


	// ------------------------ Node definition ------------------------
	// forward declaration for "friend"
	class EvidOcTree;

	class EvidOcTreeNode : public OcTreeNode {
	public:
		friend class EvidOcTree;

		EvidOcTreeNode() : updatedTime(0, 0) {
			massPtr = new EvidMass();
		}

		~EvidOcTreeNode() {
			delete massPtr;
		}

		void copyData(const EvidOcTreeNode& from) {
			// copy data payload
			OcTreeNode::copyData(from);
			// copy mass
			copyMass(from.getMassPtr());
			// copy updated time
			updatedTime = from.updatedTime;
		}
		
		/// Mass handler
		inline EvidMass* getMassPtr() const {return massPtr;}

		inline void copyMass(const EvidMass* srcMassPtr) {
			massPtr->setValue(srcMassPtr->c(), srcMassPtr->f(), srcMassPtr->o(), srcMassPtr->i());
		}

		inline void decayMass(float alp) {
			massPtr->setValue(massPtr->c() * alp,
												massPtr->f() * alp,
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

	protected:
		EvidMass* massPtr;
		ros::Time updatedTime;
	};
	 

	// ------------------------ Tree definition ------------------------
	class EvidOcTree : public OccupancyOcTreeBase<EvidOcTreeNode> {
	public:
		/// Default constructor: set resolution of leaves
		EvidOcTree(double resolution);
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
		
		void setIncomingTime(const ros::Time& _time) {incomingTime = _time;}

		/// conflict_cells_center handle
		std::vector<point3d>& getConflictCells() {return conflict_cells_center;} 
		void clearConfictCells() {conflict_cells_center.clear();}

	protected:
		// Evidential Fusion constants
		const float tau = 1.3f;  // time constant
		const float lambda_occupied = 0.7f;
		const float lambda_free = 0.7f;
		const float conflict_thres = 0.5f;  // 0.35f
		std::vector<point3d> conflict_cells_center;

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
