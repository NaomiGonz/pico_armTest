#ifndef RRT_BASE_HPP
#define RRT_BASE_HPP
#include <vector>
#include <random>
#include <cmath>
#include <stdexcept>

// Type for coordinates
using Configuration = std::vector<double>;

/*
 * Node class 
 *   - Each node represents a single point in RRT tree
 *   - Connections between nodes defined by parent/children
*/
class Node {
private:
    static bool allClose(const Configuration& a, const Configuration& b, double rtol = 1e-5, double atol = 1e-8);

public:
    Node* parent; // Pointer to the parent node (nullptr for root)
    std::vector<Node*> children; // Pointers to child nodes
    Configuration coordinates; // Coordinates of the node
    double cost; // Cost to reach this node from the root (like length of path)

    // Constructor
    Node(Node* _parent, Configuration _coordinates, double _cost);

    // Add a child node
    void addChildren(Node* _child);

    // Gets the list of children
    const std::vector<Node*>& getChildren() const;

    // Remove a specific child 
    void removeChildren(const Configuration& _conf);

    // Sets the parent node
    void setParent(Node* _parent);

    // Sets the cost
    void setCost(double _cost);
};

class GoalBiasedGreedySteerKNeighborhoodRRTStarBase {
protected:
    std::mt19937 random_generator; // C++ random number generator
    Node* root; // Pointer to the root node
    Configuration goal_coordinates; // Goal configuration
    bool goal_found; // Flag to track if goal is reached

    // Helper function for sorting pairs (e.g., for neighborhood)
    static bool compareSecond(const std::pair<Node*, double>& a, const std::pair<Node*, double>& b);

    // Recursive helper function for finding goal path
    std::vector<std::pair<Configuration, Configuration>> look_for_goal(Node* node);

    // Recursive helper function for getting all edges
    void collect_all_edges(Node* node, std::vector<std::pair<Configuration, Configuration>>& edges);


public:
    // Constructor with seed
    GoalBiasedGreedySteerKNeighborhoodRRTStarBase(int seed);

    // Destructor 
    ~GoalBiasedGreedySteerKNeighborhoodRRTStarBase();

    // Calculates the distance between two configurations (joint angles)
    // assumes Euclidian distance in joint space
    virtual double distance(const Configuration& c1, const Configuration& c2) = 0;

    // Steers from c0 towards c, checking for collisions
    virtual Configuration steer(const Configuration& c0, const Configuration& c, double step_size) = 0; // Returns nullopt or empty vector on failure?

    // Checks if two configurations are close enough
    virtual bool allclose(const Configuration& c1, const Configuration& c2) = 0;

    // Samples a configuration (with goal bias)
    virtual Configuration sample(double p) = 0;

    // Checks if a configuration is valid (collision-free)
    virtual bool valid(const Configuration& c) = 0;

    // Checks if the path between two configurations is collision-free
    virtual bool collision_free(const Configuration& c1, const Configuration& c2, double step_size) = 0;

    // --- RRT* Core Functions ---

    // Finds the k nearest neighbors to configuration c
    std::vector<Node*> neighborhood(const Configuration& c, int k);

    // Initializes/resets the RRT 
    // Creates the root node at the initial configuration (coordinate) c_init and stores the c_goal
    void init_rrt(const Configuration& c_init, const Configuration& c_goal);



    // Adds a node to the tree using RRT* logic
    // 1. Call sample(p) to get a random configuration (new_c), potentially biased towards the goal
    // 2. Call neighborhood(new_c, 1) to find the absolute nearest node (x_nearest) already in the tree
    // 3. Call steer(...) to try and move from x_nearest towards new_c by step_size, getting new_node_c
    // 4. Check collision_free between x_nearest and new_node_c
    // 5. (OPTINAL) RRT* Specific - Finding Best Parent: Calls neighborhood(new_node_c, k) to find k nearby nodes. It iterates through these near_ones to find the one (x_min) that can connect to new_node_c collision-free with the lowest total cost from the root (minimize path length)
    // 6. Adds the new_node to the tree as a child of x_min
    // 7. (OPTINAL) RRT* Specfic -  Iterates through the near_ones again. If connecting new_node to a near node results in a shorter path for that near node than its current path, it rewires the tree: disconnects near from its old parent and connects it as a child of new_node
    Configuration add_node(double p, int k, double step_size);

    
    // Gets the cost of a node (cost from root)
    double Cost(Node* node);
    
    // Gets the path from the root to the goal
    std::vector<std::pair<Configuration, Configuration>> get_path_to_goal();

    // Gets the points from the root to the goal
    std::vector<Configuration> get_points_to_goal();

    // Checks if the goal is currently reachable
    bool is_goal_reachable();

    // Simplifies a given path by trying to connect non-adjacent nodes directly if the connection is collision_free, removing intermediate nodes to shorten the path
    std::vector<std::pair<Configuration, Configuration>> simplify_path(std::vector<std::pair<Configuration, Configuration>> path, double step_size);

    // Gets all edges in the tree
    std::vector<std::pair<Configuration, Configuration>> get_all_edges();

    // Calculates the total length of the path returned by get_path_to_goal() using the distance function
    double get_goal_cost();

};

#endif //RRT_BASE_HPP