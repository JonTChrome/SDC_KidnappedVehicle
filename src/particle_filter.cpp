/*
 * particle_filter.cpp
 *
 *  Created on: Dec 12, 2016
 *      Author: Tiffany Huang
 */

#include <random>
#include <algorithm>
#include <iostream>
#include <numeric>
#include <math.h>
#include <iostream>
#include <sstream>
#include <string>
#include <iterator>

#include "particle_filter.h"

using namespace std;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
    // TODO: Set the number of particles. Initialize all particles to first position (based on estimates of
    //   x, y, theta and their uncertainties from GPS) and all weights to 1.
    // Add random Gaussian noise to each particle.
    // NOTE: Consult particle_filter.h for more information about this method (and others in this file).
    
    // Initializing the number of particles
    
    default_random_engine gen;
    num_particles = 6;
    
    double std_x = std[0];
    double std_y = std[1];
    double std_theta = std[2];
    
    normal_distribution<double> dist_x(x, std_x);
    normal_distribution<double> dist_y(y, std_y);
    normal_distribution<double> dist_theta(theta, std_theta);
    
    particles.resize(num_particles);
    weights.resize(num_particles, 1.0);
    
    for (int i = 0; i < num_particles; i++) {
        particles[i].id = i;
        particles[i].x = dist_x(gen);
        particles[i].y = dist_y(gen);
        particles[i].theta = dist_theta(gen);
        particles[i].weight = 1.0;
    }
    
    // we are initialized.
    is_initialized = true;
}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
    // TODO: Add measurements to each particle and add random Gaussian noise.
    // NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
    //  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
    //  http://www.cplusplus.com/reference/random/default_random_engine/
    
    default_random_engine gen;

    double std_x = std_pos[0];
    double std_y = std_pos[1];
    double std_theta = std_pos[2];
    
    normal_distribution<double> dist_x(0, std_x);
    normal_distribution<double> dist_y(0, std_y);
    normal_distribution<double> dist_theta(0, std_theta);
    
    //
    for (int i = 0; i < num_particles; i++) {
        
        double theta = particles[i].theta;
        
        if (fabs(yaw_rate) < 0.00001) {
            particles[i].x += velocity * delta_t * cos(theta);
            particles[i].y += velocity * delta_t * sin(theta);
        } else {
            particles[i].x += velocity / yaw_rate * (sin(theta + yaw_rate * delta_t) - sin(theta));
            particles[i].y += velocity / yaw_rate * (cos(theta) - cos(theta + yaw_rate * delta_t));
            particles[i].theta += yaw_rate * delta_t;
        }
        
        //noise
        particles[i].x += dist_x(gen);
        particles[i].y += dist_y(gen);
        particles[i].theta += dist_theta(gen);
    }
}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
    // TODO: Find the predicted measurement that is closest to each observed measurement and assign the
    //   observed measurement to this particular landmark.
    // NOTE: this method will NOT be called by the grading code. But you will probably find it useful to
    //   implement this method and use it as a helper during the updateWeights phase.

    for (unsigned int i = 0; i < observations.size(); i++) {
        //largest distance
        double minDistance = numeric_limits<double>::max();
        
        int associated_id = -1;
        
        for (unsigned j = 0; j < predicted.size(); j++) {
            
            double distance = dist(observations[i].x, observations[i].y, predicted[j].x, predicted[j].y);
            
            if (distance < minDistance) {
                minDistance = distance;
                associated_id = predicted[j].id;
            }
        }
        
        observations[i].id = associated_id;
    }
}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[],
                                   const std::vector<LandmarkObs> &observations, const Map &map_landmarks) {
    // TODO: Update the weights of each particle using a mult-variate Gaussian distribution. You can read
    //   more about this distribution here: https://en.wikipedia.org/wiki/Multivariate_normal_distribution
    // NOTE: The observations are given in the VEHICLE'S coordinate system. Your particles are located
    //   according to the MAP'S coordinate system. You will need to transform between the two systems.
    //   Keep in mind that this transformation requires both rotation AND translation (but no scaling).
    //   The following is a good resource for the theory:
    //   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
    //   and the following is a good resource for the actual equation to implement (look at equation
    //   3.33
    //   http://planning.cs.uiuc.edu/node99.html
    double range = std_landmark[0];
    double bearing = std_landmark[1];
    
    for (int i = 0; i < num_particles; i++) {
        
        double x = particles[i].x;
        double y = particles[i].y;
        double theta = particles[i].theta;
        
        double sensor_range_2 = pow(sensor_range, 2);
        
        vector<LandmarkObs> landmarks;
        
        for(int j = 0; j < map_landmarks.landmark_list.size(); j++) {
            float lm_x = map_landmarks.landmark_list[j].x_f;
            float lm_y = map_landmarks.landmark_list[j].y_f;
            int lm_id = map_landmarks.landmark_list[j].id_i;
            
            if (pow((x - lm_x), 2) + pow((y - lm_y), 2) <= sensor_range_2) {
                landmarks.push_back(LandmarkObs{lm_id, lm_x, lm_y});
            }
        }
        
        vector<LandmarkObs> transformed;
        for(int j = 0; j < observations.size(); j++) {
            double tr_x = cos(theta)*observations[j].x - sin(theta)*observations[j].y + x;
            double tr_y = sin(theta)*observations[j].x + cos(theta)*observations[j].y + y;
            transformed.push_back(LandmarkObs{observations[j].id, tr_x, tr_y});
        }
        
        dataAssociation(landmarks, transformed);
        
        particles[i].weight = 1.0;
        weights[i] = 1.0;
        
        for(unsigned int j = 0; j < transformed.size(); j++) {
            
            double obs_x= transformed[j].x;
            double obs_y = transformed[j].y;
            int obs_id = transformed[j].id;
            
            double lm_x, lm_y;
            
            for(int k = 0; k < landmarks.size(); k++) {
                if (landmarks[k].id == obs_id) {
                    lm_x = landmarks[k].x;
                    lm_y = landmarks[k].y;
                    break;
                }
            }
            
            double dX = obs_x - lm_x;
            double dY = obs_y - lm_y;
            
            double sig_x = std_landmark[0];
            double sig_y = std_landmark[1];
            double sig_y_2 = pow(sig_x, 2);
            double sig_x_2 = pow(sig_y, 2);
            double normalizer = (1.0/(2.0 * M_PI * sig_x * sig_y));
            
            double weight = normalizer * exp(-(pow(dX, 2)/(2.0 * sig_x_2) + (pow(dY, 2)/(2.0 * sig_y_2))));
            if (weight == 0) {
                particles[i].weight *= 0.00001;
                weights[i] *= 0.000001;
            } else {
                particles[i].weight *= weight;
                weights[i] *= weight;
            }
        }
    }
}

void ParticleFilter::resample() {
    // TODO: Resample particles with replacement with probability proportional to their weight.
    // NOTE: You may find std::discrete_distribution helpful here.
    //   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
    
    default_random_engine gen;
    
    double maxWeight = numeric_limits<double>::min();
    for(int i = 0; i < num_particles; i++) {
        if (particles[i].weight > maxWeight) {
            maxWeight = particles[i].weight;
        }
    }
    
    uniform_real_distribution<double> distDouble(0.0, maxWeight);
    uniform_int_distribution<int> distInt(0, num_particles - 1);
    
    int index = distInt(gen);
    double beta = 0.0;
    
    vector<Particle> resampled;
    for(int i = 0; i < num_particles; i++) {
        beta += distDouble(gen) * 2.0;
        while(beta > weights[index]) {
            beta -= weights[index];
            index = (index + 1) % num_particles;
        }
        resampled.push_back(particles[index]);
    }
    
    particles = resampled;
}

Particle ParticleFilter::SetAssociations(Particle& particle, const std::vector<int>& associations,
                                         const std::vector<double>& sense_x, const std::vector<double>& sense_y)
{
    //particle: the particle to assign each listed association, and association's (x,y) world coordinates mapping to
    // associations: The landmark id that goes along with each listed association
    // sense_x: the associations x mapping already converted to world coordinates
    // sense_y: the associations y mapping already converted to world coordinates
    
    //Clear the previous associations
    particle.associations.clear();
    particle.sense_x.clear();
    particle.sense_y.clear();
    
    particle.associations= associations;
    particle.sense_x = sense_x;
    particle.sense_y = sense_y;
    
    return particle;
}

string ParticleFilter::getAssociations(Particle best)
{
    vector<int> v = best.associations;
    stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<int>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseX(Particle best)
{
    vector<double> v = best.sense_x;
    stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseY(Particle best)
{
    vector<double> v = best.sense_y;
    stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
