#include "rrt_base.hpp"
#include <vector>
#include <cmath>
#include <algorithm> 
#include <limits>    
#include <stdexcept> 
#include <queue>     
#include <map>  


Node::Node(Node* _parent, Configuration _coordinates, double _cost)
    : parent(_parent),
      coordinates(std::move(_coordinates)), 
      cost(_cost)
{
    
}

void Node::addChildren(Node* _child) {
    children.push_back(_child);
}

const std::vector<Node*>& Node::getChildren() const {
    return children;
}

bool Node::allClose(const Configuration& a, const Configuration& b, double rtol, double atol)
{
    if (a.size() != b.size()) return false;

    for (std::size_t i = 0; i < a.size(); ++i) {
        if (std::fabs(a[i] - b[i]) > (atol + rtol * std::fabs(b[i])))
            return false;
    }
    return true;
}

void Node::removeChildren(const Configuration& _conf)
{
    children.erase(
        std::remove_if(children.begin(), children.end(),[&](Node* child) {
            return allClose(child->coordinates, _conf);
        }),
        children.end());
}

void Node::setParent(Node* _parent) {
    parent = _parent;
}

void Node::setCost(double _cost) {
    cost = _cost;
}

GoalBiasedGreedySteerKNeighborhoodRRTStarBase::GoalBiasedGreedySteerKNeighborhoodRRTStarBase(int seed)
    : root(nullptr), 
      goal_found(false) 
{
    random_generator.seed(seed); 
}

GoalBiasedGreedySteerKNeighborhoodRRTStarBase::~GoalBiasedGreedySteerKNeighborhoodRRTStarBase() {
    // Need to delete all nodes in the tree to prevent memory leaks
    if (root) {
        std::vector<Node*> nodes_to_delete;
        std::queue<Node*> queue;
        queue.push(root);

        while (!queue.empty()) {
            Node* current = queue.front();
            queue.pop();
            nodes_to_delete.push_back(current); // Add node to delete list

            for (Node* child : current->getChildren()) {
                queue.push(child); // Add children to queue
            }
        }

        // Delete all collected nodes
        for (Node* node : nodes_to_delete) {
            delete node;
        }
        root = nullptr; // Ensure root is null after deletion
    }
}

// Helper function for sorting pairs by the second element (distance)
bool GoalBiasedGreedySteerKNeighborhoodRRTStarBase::compareSecond(
    const std::pair<Node*, double>& a, const std::pair<Node*, double>& b)
{
    return a.second < b.second;
}

// Finds the k nearest neighbors to configuration c using BFS
std::vector<Node*> GoalBiasedGreedySteerKNeighborhoodRRTStarBase::neighborhood(const Configuration& c, int k) {
    if (!root) {
        return std::vector<Node*>(); // Return empty vector if tree is empty
    }

    // Use a vector of pairs to store (Node*, distance)
    std::vector<std::pair<Node*, double>> neighbors;
    std::queue<Node*> queue; // Queue for Breadth-First Search (BFS)

    queue.push(root);
    neighbors.push_back({root, distance(root->coordinates, c)}); // Add root initially

    while (!queue.empty()) {
        Node* current_node = queue.front();
        queue.pop();

        // Process children
        for (Node* child : current_node->getChildren()) {
            if (!child) continue; // Safety check

            double dist_to_c = distance(child->coordinates, c);
            queue.push(child); // Add child to queue for further traversal

            // Maintain a sorted list of the k nearest neighbors found so far
            if (neighbors.size() < k) {
                // If we haven't found k neighbors yet, just add and sort
                neighbors.push_back({child, dist_to_c});
                std::sort(neighbors.begin(), neighbors.end(), compareSecond);
            } else if (dist_to_c < neighbors.back().second) {
                // If the new node is closer than the furthest neighbor currently stored
                neighbors.pop_back(); // Remove the furthest
                neighbors.push_back({child, dist_to_c}); // Add the new one
                std::sort(neighbors.begin(), neighbors.end(), compareSecond); // Re-sort
            }
        }
    }

    // Extract just the Node pointers from the sorted pairs
    std::vector<Node*> result_nodes;
    result_nodes.reserve(neighbors.size());
    for (const auto& pair : neighbors) {
        result_nodes.push_back(pair.first);
    }

    return result_nodes;
}


// Initializes/resets the RRT
void GoalBiasedGreedySteerKNeighborhoodRRTStarBase::init_rrt(const Configuration& c_init, const Configuration& c_goal) {
    // Clean up existing tree if it exists
    if (root) {
        // Use the destructor logic to delete all nodes
        std::vector<Node*> nodes_to_delete;
        std::queue<Node*> queue;
        queue.push(root);
        while (!queue.empty()) {
             Node* current = queue.front(); queue.pop();
             nodes_to_delete.push_back(current);
             for (Node* child : current->getChildren()) { queue.push(child); }
        }
        for (Node* node : nodes_to_delete) { delete node; }
    }

    // Create the new root node
    root = new Node(nullptr, c_init, 0.0); // Parent is null, cost is 0

    // Store goal coordinates and reset goal flag
    goal_coordinates = c_goal;
    goal_found = false; // Reset goal found status

     // Check if init configuration itself is close to goal
     if (root && allclose(root->coordinates, goal_coordinates)) {
         goal_found = true;
     }
}

// Gets the cost of a node from root
double GoalBiasedGreedySteerKNeighborhoodRRTStarBase::Cost(Node* node) {
    if (!node) {
        return std::numeric_limits<double>::infinity();
    }
    return node->cost;
}

// Checks if the goal is currently reachable
bool GoalBiasedGreedySteerKNeighborhoodRRTStarBase::is_goal_reachable() {
    if(goal_found) return true;

    // Verify by checking path existence 
    std::vector<std::pair<Configuration, Configuration>> path = get_path_to_goal();
    if (!path.empty()) {
         goal_found = true; 
         return true;
    }
    return false;

}

// Adds a node to the tree using RRT* logic
Configuration GoalBiasedGreedySteerKNeighborhoodRRTStarBase::add_node(double p, int k, double step_size) {
    const int MAX_SAMPLE_ATTEMPTS = 20; // Limit attempts to find a steerable point
    Configuration new_node_c;         // Will hold the config of the node to add
    Node* x_nearest = nullptr;        // Nearest node found
    bool steer_success = false;

    // Loop to find a valid point using sampling and steering
    for (int attempt = 0; attempt < MAX_SAMPLE_ATTEMPTS; ++attempt) {
        // Sample a configuration
        Configuration new_c = sample(p);
        if (new_c.empty()) {
            continue; 
        }

        // Find the nearest node in the tree as potential parent 
        std::vector<Node*> nearest_nodes = neighborhood(new_c, 1);
        if (nearest_nodes.empty()) {
            // Tree empty, cant add node
             return Configuration(); 
        }
        x_nearest = nearest_nodes[0]; 

        // Steer from nearest node towards the sampled configuration
        new_node_c = steer(x_nearest->coordinates, new_c, step_size);

        // Don't add if the steered point is identical to the nearest node
        if (!new_node_c.empty()) {
            if (!allclose(new_node_c, x_nearest->coordinates)) {
                 steer_success = true; 
                 break; 
            }
        }
        // If steering failed loop continues to try again
    } 

    // Check if we successfully found a node 
    if (!steer_success || x_nearest == nullptr) {
        return Configuration(); 
    }

    // --- RRT* Specific Logic  ---

    // Find k nearest neighbors to the new potential node 
    std::vector<Node*> near_neighbors = neighborhood(new_node_c, k);

    // Find the best parent which minimize cost
    Node* x_min = x_nearest;
    double c_min = Cost(x_nearest) + distance(x_nearest->coordinates, new_node_c);

    // Check if path (x_nearest -> new_node_c) is valid
    if (!collision_free(x_nearest->coordinates, new_node_c, step_size)) {
         c_min = std::numeric_limits<double>::infinity(); 
         x_min = nullptr; 
    }

    // Iterate through neighbors to find the best valid parent which is collision free
    for (Node* near_node : near_neighbors) {
        if (collision_free(near_node->coordinates, new_node_c, step_size)) {
            double cost_via_near = Cost(near_node) + distance(near_node->coordinates, new_node_c);
            if (cost_via_near < c_min) {
                c_min = cost_via_near; 
                x_min = near_node;     
            }
        }
    }

    // If a valid parent not found return empyt vector
    if (x_min == nullptr) {
        return Configuration(); 
    }

    // Add the new node to the tree
    Node* new_node = new Node(x_min, new_node_c, c_min);
    x_min->addChildren(new_node); 

    // Check if new node offers a better path for its neighbors
    for (Node* near_node : near_neighbors) {
        if (near_node == x_min) {
            continue;
        }

        // Calculate cost to reach near_node using the new_node
        double cost_via_new = Cost(new_node) + distance(new_node_c, near_node->coordinates);

        // Check if shorter path using new_node and collsion free 
        if (cost_via_new < Cost(near_node) && collision_free(new_node_c, near_node->coordinates, step_size)) {

        Node* old_parent = near_node->parent;
        if (old_parent) {
            // Detach near_node from its previous parent
            old_parent->removeChildren(near_node->coordinates);
        }

        near_node->setParent(new_node);
        near_node->setCost(cost_via_new);
        new_node->addChildren(near_node);
        }
    }

    // Check if the newly added node reached the goal
    if (allclose(new_node->coordinates, goal_coordinates)) {
        goal_found = true; // Set the flag
    }

    // Return the configuration of the newly added node
    return new_node_c;
}


// Simplifies a given path
std::vector<std::pair<Configuration, Configuration>> GoalBiasedGreedySteerKNeighborhoodRRTStarBase::simplify_path(
    std::vector<std::pair<Configuration, Configuration>> path, double step_size)
{
    if (path.size() <= 1) {
        return path; // Cannot simplify path with 0 or 1 segments
    }

    int i = 0;
    while (i < path.size()) {
        // Start point of the current segment (edge i)
        const Configuration& start_node_config = path[i].first;
        int best_reachable_j = -1; // Index of the *end* node of the furthest reachable segment

        // Try to connect start_node_config to the end nodes of subsequent segments
        // Check from furthest possible end node (path.back().second) backwards
        for (int j = path.size() - 1; j > i; --j) {
             const Configuration& potential_end_node_config = path[j].second;

             if (collision_free(start_node_config, potential_end_node_config, step_size)) {
                 // Found a direct connection from start_node_config to the end of segment j
                 best_reachable_j = j;
                 break; // Stop searching backwards once the furthest clear connection is found
             }
        }

        // If a shortcut was found (best_reachable_j > i)
        if (best_reachable_j > i) {
            // Modify the current segment 'i' to connect directly to the end of segment 'best_reachable_j'
            path[i].second = path[best_reachable_j].second;

            // Remove the intermediate segments (from i+1 up to and including best_reachable_j)
            path.erase(path.begin() + i + 1, path.begin() + best_reachable_j + 1);

            // Don't increment 'i', re-evaluate shortcuts from the modified segment 'i'
        } else {
            // No shortcut found starting from segment 'i', move to the next segment
            i++;
        }
    }

    return path;
}

std::vector<std::pair<Configuration, Configuration>> 
GoalBiasedGreedySteerKNeighborhoodRRTStarBase::get_path_to_goal(){
    return look_for_goal(root);
}

std::vector<std::pair<Configuration, Configuration>> 
GoalBiasedGreedySteerKNeighborhoodRRTStarBase::look_for_goal(Node* node) {
    std::vector<Node*> children = node->getChildren();

    // case: goal is found directly under this node
    for(Node* child : children){
        if(allclose(child->coordinates, goal_coordinates)){
            std::vector<std::pair<Configuration, Configuration>> result;
            result.emplace_back(node->coordinates, child->coordinates);
            return result;
        }

        // case: goal not directly under this node, recurse through the subtree
        std::vector<std::pair<Configuration, Configuration>> subpath = look_for_goal(child);
        if(!subpath.empty()){
            subpath.insert(subpath.begin(), std::make_pair(node->coordinates, child->coordinates));
            return subpath;
        }
    }

    // case: goal not found in this branch
    return std::vector<std::pair<Configuration, Configuration>>();
}

std::vector<Configuration> get_points_to_goal(){
    std::vector<std::pair<Configuration, Configuration>> path = get_path_to_goal();
    std::vector<Configuration> result;

    // Loop through and only append first config
    for (int i = 0; i < path.size(); ++i){
        result.push_back(path[i].first);
    }

    return result;
}

std::vector<std::pair<Configuration, Configuration>> 
GoalBiasedGreedySteerKNeighborhoodRRTStarBase::get_all_edges(){
    std::vector<std::pair<Configuration, Configuration>> edges;

    // if root is null, return empty edge list
    if(!root) return edges;

    // use a queue for traversal and initialize with root node
    std::vector<Node*> queue;
    queue.push_back(root);

    // loop until all nodes have been visited
    while(!queue.empty()){
        //take last node from queue
        Node* current = queue.back();
        queue.pop_back(); // remove

        // get all children of the current node
        const std::vector<Node*>& children = current->getChildren();
        
        // for each child
        for(Node* child : children){
            // add the edge from teh current to child to the edge list
            edges.emplace_back(current->coordinates, child->coordinates);
            //add the child to the queue to visit its children later
            queue.push_back(child);
        }
    }

    //return list of all edges found in the tree
    return edges;
}

double GoalBiasedGreedySteerKNeighborhoodRRTStarBase::get_goal_cost(){
    //get current path from root to goal
    std::vector<std::pair<Configuration, Configuration>> path = get_path_to_goal();

    //if goal is not reachable, return infinity
    if(path.empty()){
        return std::numeric_limits<double>::infinity();
    }

    //initialize cost to 0
    double total_cost = 0.0;

    for(const auto& edge : path){
        const Configuration& from = edge.first;
        const Configuration& to = edge.second;
        total_cost += distance(from, to); // use distance() from subclass
    }

    //return cost of the path
    return total_cost;
}