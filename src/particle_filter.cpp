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

        /********************************************************************************************/ 
        //
        // Entered code 08-20-2017 cam
        //
        /********************************************************************************************/

        default_random_engine gen;

	// This line creates a normal (Gaussian) distribution for x
	normal_distribution<double> dist_x  (x,     std[0]);
	normal_distribution<double> dist_y  (y,     std[1]);
	normal_distribution<double> dist_psi(theta, std[2]);

        // number of particles
        num_particles = 100;
        
        // set all vector sizes
        weights.resize(num_particles,1.0);
        particles.resize(num_particles);

        // create particles 	
	for (int i = 0; i < num_particles ; ++i) {

	  particles[i].id      = i;
	  particles[i].x       = dist_x(gen);
	  particles[i].y       = dist_y(gen);
	  particles[i].theta   = dist_psi(gen);
	  particles[i].weight  = 1;
        
	}
        
        // set initialized
        is_initialized = true;

        cout << "called init" << " using number of particles = " << num_particles << endl << endl;
}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/

        /********************************************************************************************/ 
        //
        // Entered code 08-20-2017 cam
        //
        /********************************************************************************************/

        default_random_engine gen;

        const float ZERO_LIMIT = 0.0001;

        // predicted values
        double x;
        double y;
        double r;


        for (int i; i < num_particles; ++i) { 

          // predict next x, y, theta
          if (abs(yaw_rate) < ZERO_LIMIT) {

            x = particles[i].x + velocity * delta_t * cos(particles[i].theta);
  	    y = particles[i].y + velocity * delta_t * sin(particles[i].theta);
            r = particles[i].theta;
          }
          else {    

            x = particles[i].x + velocity / yaw_rate * ( sin(particles[i].theta + yaw_rate * delta_t) - sin(particles[i].theta));
  	    y = particles[i].y + velocity / yaw_rate * ( cos(particles[i].theta) - cos(particles[i].theta + yaw_rate * delta_t));
            r = particles[i].theta + yaw_rate * delta_t;
          } 

	// This line creates a normal (Gaussian) distribution for adding noise 
	normal_distribution<double> dist_x  (x,     std_pos[0]);
	normal_distribution<double> dist_y  (y,     std_pos[1]);
	normal_distribution<double> dist_psi(r,     std_pos[2]);
       
	particles[i].x       = dist_x(gen);
	particles[i].y       = dist_y(gen);
	particles[i].theta   = dist_psi(gen);
 
     }
       
     // cout << "called predict" << endl;

}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.

        /********************************************************************************************/ 
        //
        // Entered code 08-21-2017 cam
        //
        /********************************************************************************************/

        // predicted means:   predicted measurements for between one  particle and all map landmarks within the sensor ranges        
        // observation means: all map landmarks measurements from the lidar -> 
        //                    perform nearest data association, assign sensor assiciation          

        const double HIGH_NUM = 1000000.0;
 
        // check for each landmark measured by lidar
        for (int i=0; i<observations.size();++i) {
         
           // run through all landmarks in range and get smallest distance
           double min_dis = HIGH_NUM; 
           for (int u=0;u<predicted.size();++u) {

              double dis = dist(observations[i].x,observations[i].y,predicted[u].x,predicted[u].y);
              
              if (dis <  min_dis) {
                min_dis = dis;
                observations[i].id  = u;  
              }
           }
           // check if min_dis has been found
           if (min_dis == HIGH_NUM) cout << " Error: no min distance found for observation =" << i << endl;
             
        }       
}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
		std::vector<LandmarkObs> observations, Map map_landmarks) {
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

        /********************************************************************************************/ 
        //
        // Entered code 08-20-2017 cam
        //
        /********************************************************************************************/

        // 1. predict landmark measurements within the sensorrange for each particle
        // 2. run the dataAssioaction func. to get sensro measurements to each landmark
        // 3. calc. gaussian
        // 4. normalize weights - done in resample method using   

        // for each particle
        for (int i=0;i<num_particles;++i) {

          ///////////////////////////////////////////////////////////////////////////////
          // filtered landmarks
          ///////////////////////////////////////////////////////////////////////////////
          std::vector<LandmarkObs> landmarks_map;
          LandmarkObs lm;

          //  filter landsmarks which are in sensor range
          for (int l=0;l<map_landmarks.landmark_list.size();++l) {

             double x_f = map_landmarks.landmark_list[l].x_f;
             double y_f = map_landmarks.landmark_list[l].y_f;
         
             // check if in sensor range             
             if (dist(particles[i].x,particles[i].y,x_f,y_f) < sensor_range) {

               lm.id = map_landmarks.landmark_list[l].id_i;
               lm.x  = x_f;
               lm.y  = y_f;

               landmarks_map.push_back(lm);
             }           
          }

          ///////////////////////////////////////////////////////////////////////////////
          // transform from vehcile to map coord.
          ///////////////////////////////////////////////////////////////////////////////
          std::vector<LandmarkObs> observations_map;
          observations_map.resize(observations.size());
        
          for (int k=0;k<observations.size();++k) {

             // counterwise rotation  
             observations_map[k].id  =  observations[k].id;         
             observations_map[k].x   =  particles[i].x + cos(particles[i].theta) * observations[k].x - sin(particles[i].theta) * observations[k].y;          
             observations_map[k].y   =  particles[i].y + sin(particles[i].theta) * observations[k].x + cos(particles[i].theta) * observations[k].y; 
          }

          ///////////////////////////////////////////////////////////////////////////////
          // run data association
          ///////////////////////////////////////////////////////////////////////////////
          dataAssociation(landmarks_map,observations_map); 


          ///////////////////////////////////////////////////////////////////////////////
          // calculate weights 
          ///////////////////////////////////////////////////////////////////////////////
          double gauss_norm = (1/(2 * M_PI * std_landmark[0] * std_landmark[1]));
          double weight     = 1.0;

          for (int w=0;w<observations_map.size();++w) {

             int    index    =   observations_map[w].id; // index may be different to org. landmark id because of filter landmarks map

             double exponent = ((observations_map[w].x - landmarks_map[index].x) * (observations_map[w].x - landmarks_map[index].x)) / (2 * std_landmark[0] * std_landmark[0]) + 
                               ((observations_map[w].y - landmarks_map[index].y) * (observations_map[w].y - landmarks_map[index].y)) / (2 * std_landmark[1] * std_landmark[1]);

             weight *= gauss_norm * exp(-exponent); 

          }

          // assign to particle
          particles[i].weight = weight;
          weights[i]          = weight;

          // empty vectors
          observations_map.clear();
          landmarks_map.clear();

      }
      
      // cout << "called update" << endl;
}


void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
 
        ///////////////////////////////////////////////////////////////////////////////
        // distribute according to above 
        ///////////////////////////////////////////////////////////////////////////////
        default_random_engine gen;
        discrete_distribution<int> distribution (weights.begin(),weights.end());
 
        // setup resample vector      
        std::vector<Particle> next_particles;

         // resample 
        for (int i=0;i<num_particles;++i) {
           int index = distribution(gen);
           next_particles.push_back(particles[index]);
        }

        // set particles
        particles = next_particles;    

        //cout << "called resample" << endl << endl;
}

Particle ParticleFilter::SetAssociations(Particle particle, std::vector<int> associations, std::vector<double> sense_x, std::vector<double> sense_y)
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
