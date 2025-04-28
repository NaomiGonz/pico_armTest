#ifndef RRT_STAR_3D_H
#define RRT_STAR_3D_H

#include "rrt_base.hpp" 
#include <vector>
#include <cmath>
#include <stdexcept> 

// Structure for sphere obstacles
struct SphereObstacle {
    double x, y, z, radius;
};

// Should we make a structure that holds arm dimenstions as a rectangle?
// or we treat it as a line segment?

class RRTStar3D : public GoalBiasedGreedySteerKNeighborhoodRRTStarBase {
private:
    std::vector<SphereObstacle> obstacles;  
    double xyz_min, xyz_max;
    int dof;
    std::vector<std::pair<double, double>> rad_limits;                           

    // Takes a joint configuration and returns all joint points in Cartesian space? or just final point?
    std::vector<Configuration> forward_kinematics(const Configuration& joint_angles);

    // Check if robot arm collidied with spheres in cartesian (after forward kinematics)
    bool check_robot_collision(const std::vector<Configuration>& robot_cartesian_points);

    // similar to the 2D intersect but in 3D if each arm is a "line segment"
    bool segment_sphere_intersect(const Configuration& p1, const Configuration& p2, const SphereObstacle& obs);

public:

    RRTStar3D(int seed,
              const std::vector<SphereObstacle>& _obstacles,
              const std::vector<std::pair<double, double>>& _rad_limits
             );

    // --- Implementations of Virtual Functions from Base Class ---

    // Gidon
    // Calculates distance between two joint configurations (6D vectors)
    double distance(const Configuration& c1, const Configuration& c2) override;

    // Naomi
    // Steers from joint config c0 towards c by a joint-space step_size
    // checks collision free using FK at each step taken
    Configuration steer(const Configuration& c0, const Configuration& c, double step_size) override;

    // Gidon 
    // use joint-space distance and check closeness
    bool allclose(const Configuration& c1, const Configuration& c2) override;


    // Gidon 
    // Sample random 6-DOF joint configuration within limits
    Configuration sample(double p) override;

    // Naomi
    // Check if valid configuration or if it hit obstacles 
    bool valid(const Configuration& c) override;

    // Naomi
    // check if path in joint space is collision free
    bool collision_free(const Configuration& c1, const Configuration& c2, double step_size) override;
};

#endif