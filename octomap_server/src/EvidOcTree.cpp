#include <octomap_server/EvidOcTree.h>
#include <cmath>

namespace octomap {

// ------------------------ Node impl ------------------------
void EvidOcTreeNode::updateOccupancyChildren() {
	// Compute children's average mass
	float f_(0.0f), o_(0.0f), i_(0.0f);
	int num_child(0);
	
	if(children != NULL) {
		for(int i=0; i<8; i++) {
			EvidOcTreeNode* child = static_cast<EvidOcTreeNode*>(children[i]);
			if(child != NULL) {
				f_ += child->massPtr->f();
				o_ += child->massPtr->o();
				i_ += child->massPtr->i();
			}
		}
	}

	if(num_child > 0) {
		f_ /= (float) num_child;
		o_ /= (float) num_child;
		i_ /= (float) num_child;
	}

	// set node's evid mass
	massPtr->setValue(f_, o_, i_);

	// convert evid mass to log odds
	setLogOdds(evidMassToLogOdds());
}



// ------------------------ Tree impl ------------------------
EvidOcTree::StaticMemberInitializer EvidOcTree::evidOcTreeMemberInit;

EvidOcTree::EvidOcTree(double in_resolution) :
	OccupancyOcTreeBase<EvidOcTreeNode>(in_resolution),
	basic_belief_occupied(0.0f, lambda_occupied, 1.0f-lambda_occupied),
	basic_belief_free(lambda_free, 0.0f, 1.0f-lambda_free)
	{
		evidOcTreeMemberInit.ensureLinking();
	};


EvidOcTreeNode* EvidOcTree::updateNode(const OcTreeKey& key, bool occupied, bool lazy_eval) {
	if(occupied) {
		return updateNode(key, basic_belief_occupied, lazy_eval);
	} else {
		return updateNode(key, basic_belief_free, lazy_eval);
	}
	
	// Rewire the flow of node update
	// return updateNode(key, basicBeliefAssign, lazy_eval);
}


EvidOcTreeNode* EvidOcTree::updateNode(const OcTreeKey& key, const EvidMass& basicBeliefAssign, bool lazy_eval) {
	bool createdRoot = false;
	// check if tree is just initialized
	if (this->root == NULL){
		this->root = new EvidOcTreeNode();
		this->tree_size++;
		createdRoot = true;
	}

	return updateNodeRecurs(this->root, createdRoot, key, 0, basicBeliefAssign, lazy_eval);
}


EvidOcTreeNode* EvidOcTree::updateNodeRecurs(EvidOcTreeNode* node, bool node_just_created, const OcTreeKey& key,
                                    unsigned int depth, const EvidMass& basicBeliefAssign, bool lazy_eval) {
	bool created_node = false;
	assert(node);

	// follow down to last level
	if(depth < this->tree_depth) {
		unsigned int pos = computeChildIdx(key, this->tree_depth -1 - depth);
		if (!this->nodeChildExists(node, pos)) {
			// child does not exist, but maybe it's a pruned node?
			if (!this->nodeHasChildren(node) && !node_just_created ) {
				// current node does not have children AND it is not a new node
				// -> expand pruned node
				this->expandNode(node);
			} else {
				// not a pruned node, create requested child
				this->createNodeChild(node, pos);
				created_node = true;
				// init child's updatedTime
				EvidOcTreeNode* child = getNodeChild(node, pos);
				child->setUpdatedTime(incomingTime);
			}
		}

		if (lazy_eval)
			return updateNodeRecurs(this->getNodeChild(node, pos), created_node, key, depth+1, basicBeliefAssign, lazy_eval);
		else {
			EvidOcTreeNode* retval = updateNodeRecurs(this->getNodeChild(node, pos), created_node, key, depth+1, basicBeliefAssign, lazy_eval);
			// prune node if possible, otherwise set own probability
			if (this->pruneNode(node)){
			  // return pointer to current parent (pruned), the just updated node no longer exists
			  return node;  // ptr to parent
			} else {
			  node->updateOccupancyChildren();
			}
			return retval;
		}
	} 
	else { // at last level, update node, end of recursion
		upadteNodeEvidMass(node, basicBeliefAssign, key);
		return node;  // ptr to leaf
	}

}

void EvidOcTree::upadteNodeEvidMass(EvidOcTreeNode* node, const EvidMass& basicBeliefAssign, const OcTreeKey& key) {
	/// Compute decay factor (alpha)
	// time from last update of this node to the arrival of current pointcloud
	double time_elapsed = (incomingTime - node->getUpdatedTime()).toSec();
	assert(time_elapsed >= 0.0);
	// update node's updatedTime
	node->setUpdatedTime(incomingTime);

	float alpha = exp(-(float) time_elapsed / tau);
	// if (alpha == 0.0f)
		// ROS_WARN("Alpha equals to zero");

	// decay mass
	node->decayMass(alpha);

	// Conjunctive fusion (f_ means fused)
	float f_c = node->massPtr->f() * basicBeliefAssign.o() + node->massPtr->o() * basicBeliefAssign.f();
	float f_f = node->massPtr->f() * basicBeliefAssign.f() + node->massPtr->f() * basicBeliefAssign.i() + node->massPtr->i() * basicBeliefAssign.f();
	float f_o = node->massPtr->o() * basicBeliefAssign.o() + node->massPtr->o() * basicBeliefAssign.i() + node->massPtr->i() * basicBeliefAssign.o();
	float f_i = node->massPtr->i() * basicBeliefAssign.i();

	// TODO: Get 3d coord of cells having high conflict mass
	if(f_c > conflict_thres) {
		point3d _p = keyToCoord(key);
		conflict_cells_center.push_back(_p);
		// Do something
		// std::cout<<"[INFO] Detect conflict cells at (" << _p.x() << ", " << _p.y() << ", " << _p.z() <<")\n";
	} 

	// Dempster normalization
	float k = 1.0f / (1.0f - f_c);
	node->massPtr->setValue(k*f_f, k*f_o, k*f_i);
	assert(node->massPtr->isValid());

	// update node's logodds according to its evidential mass
	node->setLogOdds(node->evidMassToLogOdds());

	// clamping logodds
	if (node->getLogOdds() < this->clamping_thres_min) {
		node->setLogOdds(this->clamping_thres_min);
		return;
	}
	if (node->getLogOdds() > this->clamping_thres_max) {
		node->setLogOdds(this->clamping_thres_max);
	}
}


}// end of namespace
