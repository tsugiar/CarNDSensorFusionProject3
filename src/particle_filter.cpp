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
#include <limits>

#include "particle_filter.h"

using namespace std;


void ParticleFilter::init(double x, double y, double theta, double std[]) {
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of 
	//   x, y, theta and their uncertainties from GPS) and all weights to 1. 
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).

	std::default_random_engine generator;

	std::normal_distribution<double> x_distribution(x, std[0]);
	std::normal_distribution<double> y_distribution(y, std[1]);
	std::normal_distribution<double> theta_distribution(theta, std[2]);
	double x_sample, y_sample, theta_sample;

	// Set default size of particle filter
	num_particles = 25;
	
	// Pre-allocate memory for number of particles
	particles.reserve(num_particles);

	Particle particle_el;

	for (int i = 0; i < num_particles; i++)
	{
		is_initialized = true;

		// Reserve space for x_sense, y_sense, and data-association

		// Generate sample from random gaussian distribution based on 
		// initial position x,y, and theta
		x_sample = x_distribution(generator);
		y_sample = y_distribution(generator);
		theta_sample = theta_distribution(generator);

		// Populate particle element
		particle_el.id = i;
		particle_el.x = x_sample;
		particle_el.y = y_sample;
		particle_el.theta = theta_sample;
		particle_el.weight = 1;
	    

		particles.push_back(particle_el);

		// Reserve space for x_sense, y_sense, and data-association
		// Space is twice number of particles, to accomodate multiple detection
	
		particles.back().sense_x.reserve(num_particles * 2);
		particles.back().sense_y.reserve(num_particles * 2);
		particles.back().associations.reserve(num_particles * 2);
		particles.back().qualified_lmark.reserve(num_particles * 2);

	}


}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/

	std:default_random_engine generator;
	std::normal_distribution<double> x_dist(0, std_pos[0]);
	std::normal_distribution<double> y_dist(0, std_pos[1]);
	std::normal_distribution<double> theta_dist(0, std_pos[2]);

	double x_noise, y_noise, theta_noise;


	// Iterating over each particle
	for (int i = 0; i < particles.size(); i++)
	{
		// Generate noise for x, y, and theta
		x_noise = x_dist(generator);
		y_noise = y_dist(generator);
		theta_noise = theta_dist(generator);

		if (fabsf(yaw_rate) > 1e-5)
		{
			// Predict x
			particles[i].x = particles[i].x +  x_noise +
				velocity / yaw_rate * (sin(particles[i].theta + yaw_rate * delta_t) - sin(particles[i].theta));
			// Predict y
			particles[i].y = particles[i].y + y_noise +
				velocity / yaw_rate * (cos(particles[i].theta) - cos( particles[i].theta + yaw_rate*delta_t ));

		}
		else
		{
			particles[i].x = particles[i].x + x_noise + velocity * delta_t * cos(particles[i].theta);
			particles[i].y = particles[i].y + y_noise + velocity * delta_t * sin(particles[i].theta);

		}

		// Predict theta
		particles[i].theta = particles[i].theta + theta_noise + yaw_rate * delta_t;

	}

}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.
  

}

void ParticleFilter::AssociateObsLandmark(Particle &particle, const Map &map_landmarks, double std_landmark[])
{
	double lmark_x, lmark_y;   // Land-mark position
	double temp_dist; 
	double weight;
	double temp;
	double denom;

	struct closest_dist_s
	{ 
		// obs_index used when outer loop is land-mark
		// lmar_index used when outer loop is observation
		int obs_index;  int lmark_index;  double obs_lmark_dist; double lmark_x; double lmark_y;
	};

	closest_dist_s closest_dist;

	denom = 2 * M_PI * std_landmark[0] * std_landmark[1];

	// Iterate over qualified land-mark (stored inside particle), to determine which one 
	// is closer to sensed observation

	weight = 0.0;
    // If size of qualified land-mark less / equal than no detection, 
	// put landmark iteration on outer loop
	if (particle.qualified_lmark.size() <= particle.sense_x.size())
	{
		for (auto & qual_lmark : particle.qualified_lmark)
		{
			// qual_lmark pair variable, with first of pair is land-mark ID, second of pair is its occupation-status

			// Get x , y position of the land-mark and add-noise
			// Put qual_lmark.first - 1, because indexing on land-mark starts with index 1 instead of 0
			lmark_x = map_landmarks.landmark_list[qual_lmark.first - 1].x_f;
			lmark_y = map_landmarks.landmark_list[qual_lmark.first - 1].y_f;

			// Initialize closest_dist.second to high-number
			closest_dist.obs_lmark_dist = 1e6;

			// For each land-mark we would like to check, closest observation detection
			for (int obs_idx = 0; obs_idx < particle.sense_x.size(); ++obs_idx)
			{

				temp_dist = dist(particle.sense_x[obs_idx], particle.sense_y[obs_idx], lmark_x, lmark_y);

				// Condition to select closest observation-detection is min-distance to landmark
				// and observed-detection hasn't been associated with particular landmark
				if ((temp_dist < closest_dist.obs_lmark_dist) && (particle.associations[obs_idx] == 1e6))
				{
					closest_dist.obs_lmark_dist = temp_dist;
					closest_dist.obs_index = obs_idx;
				}

			}  // End of observation loop

			// Associate sense_x detection with landmark
			if (particle.associations[closest_dist.obs_index] == 1000000)
			{
				// Association hasn't been made (big number), then update assocations
				// and its weight vector


				particle.associations[closest_dist.obs_index] = qual_lmark.first;

				// Start updating weight,
				temp = exp(-0.5 * ((particle.sense_x[closest_dist.obs_index] - lmark_x) *
					1.0 / (std_landmark[0] * std_landmark[0])          *
					(particle.sense_x[closest_dist.obs_index] - lmark_x) +
					(particle.sense_y[closest_dist.obs_index] - lmark_y) *
					1.0 / (std_landmark[1] * std_landmark[1])          *
					(particle.sense_y[closest_dist.obs_index] - lmark_y))) / denom;

				weight = weight == 0.0 ? temp : weight * temp;

			}


		} // End of looping statement over qualified land-mark
	}
	else
	{  // If size of qualified land-mark greater than number of detection, put 
		// detection iteration on outer loop

		for (int obs_idx = 0; obs_idx < particle.sense_x.size(); ++obs_idx)
		{

			closest_dist.obs_lmark_dist = 1e6;

			for (auto & qual_lmark : particle.qualified_lmark)
			{
				lmark_x = map_landmarks.landmark_list[qual_lmark.first - 1].x_f;
				lmark_y = map_landmarks.landmark_list[qual_lmark.first - 1].y_f;


				temp_dist = dist(particle.sense_x[obs_idx], particle.sense_y[obs_idx], lmark_x, lmark_y);

				// Condition to select closest observation-detection is min-distance to landmark
				// and observed-detection hasn't been associated with particular landmark
				if ((temp_dist < closest_dist.obs_lmark_dist) && (particle.associations[obs_idx] == 1e6))
				{
					closest_dist.obs_lmark_dist = temp_dist;
					closest_dist.lmark_index = qual_lmark.first;
					closest_dist.lmark_x = lmark_x;
					closest_dist.lmark_y = lmark_y;
					closest_dist.obs_index = obs_idx;
				}


			}  // End of iterations over qualified land-mark

			   // Associate sense_x detection with landmark
			if (particle.associations[closest_dist.obs_index] == 1000000)
			{
				// Association hasn't been made (big number), then update assocations
				// and its weight vector


				particle.associations[closest_dist.obs_index] = closest_dist.lmark_index;

				// Start updating weight,
				temp = exp(-0.5 * ((particle.sense_x[closest_dist.obs_index] - closest_dist.lmark_x) *
					1.0 / (std_landmark[0] * std_landmark[0])          *
					(particle.sense_x[closest_dist.obs_index] - closest_dist.lmark_x) +
					(particle.sense_y[closest_dist.obs_index] - closest_dist.lmark_y) *
					1.0 / (std_landmark[1] * std_landmark[1])          *
					(particle.sense_y[closest_dist.obs_index] - closest_dist.lmark_y))) / denom;

				weight = weight == 0.0 ? temp : weight * temp;

			}


		}
	}  // End of else condition when size of land-mark greater than detection


	// Assign weight to particle
	particle.weight = weight;

	// Remove detections sense_x, sense_y, and its associations with land-mark whenever associations
	// is still not assigned yet (e.g. associations has high number 1e6)
	
	bool trigger_move = false;
	bool break_all_loop = false;
	int temp_next_idx = 0;   
	int temp_curr_idx = 0;

	while (temp_curr_idx < particle.sense_x.size())
	{
		if (particle.associations[temp_curr_idx] == 1e6 || trigger_move == 1)
		{
			trigger_move = true;
			
			while (particle.associations[temp_next_idx] == 1e6)
			{  // Everytime we found 1e6, increment temp_next_idx
				temp_next_idx = temp_next_idx + 1;

				if (temp_next_idx >= particle.sense_x.size())
				{
					temp_next_idx = particle.sense_x.size() - 1;
					break_all_loop = true;
					break;
				}
			}

			if (break_all_loop == true)
			{
				// also break second while loop

				break;
			}

		}

		// Update particle values
		particle.associations[temp_curr_idx] = particle.associations[temp_next_idx];
		particle.sense_x[temp_curr_idx] = particle.sense_x[temp_next_idx];
		particle.sense_y[temp_curr_idx] = particle.sense_y[temp_next_idx];

		// Increment temp_next_idx and temp_curr_idx
		temp_next_idx = temp_next_idx  < particle.sense_x.size()-1 ? temp_next_idx + 1 : temp_next_idx;
		temp_curr_idx = temp_curr_idx + 1;
	}

	// Check if temp_curr_idx is less than particle.sense_x.size(), it indicates that 
	// last few entries beginning with temp_curr_idx doesn't have land-mark association
	// Remove the associations, and sensor detection
	if (temp_curr_idx < particle.sense_x.size())
	{
		particle.associations.erase(particle.associations.begin() + temp_curr_idx, particle.associations.end());
		particle.sense_x.erase(particle.sense_x.begin() + temp_curr_idx, particle.sense_x.end());
		particle.sense_y.erase(particle.sense_y.begin() + temp_curr_idx, particle.sense_y.end());

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

	double range;
	double x_obs_world;			// Landmark observation x in map-coord system 
	double y_obs_world;			// Landmark observation y in map-coord system 
	double particle_theta;		// Particle heading angle in respect to map-coord system
	double particle_x;			// Particle x-pos in map-coord system
	double particle_y;			// Particle y-pos in map-coord system
	double cos_theta;
	double sin_theta;


	for (auto & particle_item : particles)
	{
		// Compute cos and sin theta of particle, for transformation purpose
		cos_theta = cos(particle_item.theta);
		sin_theta = sin(particle_item.theta);
		particle_x = particle_item.x;
		particle_y = particle_item.y;

		// Clear-up previous observation detection
		particle_item.sense_x.clear();   
		particle_item.sense_y.clear();
		particle_item.associations.clear();
		particle_item.qualified_lmark.clear();

		// Iterate over land-mark position, to determine if it is within sensor-range
		for (auto const lmark_item : map_landmarks.landmark_list)
		{
			range = dist(particle_item.x, particle_item.y, lmark_item.x_f, lmark_item.y_f);

			// If land-mark position within sensor-reach insert it into qualified landmark
			if (range < sensor_range )  // Add 10 m for hysteresis purpose
			{
				// Put ID of qualified land-mark and its occupied status to qualified_lmark 
				// field
				particle_item.qualified_lmark.push_back(std::make_pair(lmark_item.id_i, 0));
			}
		}


		// Iterate over each observation perceived by particles
		for (auto obs_item : observations)
		{  
			// Transformed each landmark-observation from vehicle to map coord system
			x_obs_world = cos_theta * obs_item.x - sin_theta * obs_item.y + particle_x;
			y_obs_world = sin_theta * obs_item.x + cos_theta * obs_item.y + particle_y;

			// Check range between particle to landmark-observation
			range = dist(particle_x, particle_y, x_obs_world, y_obs_world);

			// Landmark observation is valid only when its range less than / equal to sensor_range
			if (range <= sensor_range  - 5.0)  // Add 10 m for hysteresis purpose
			{
				// Process observation data-further
				particle_item.sense_x.push_back(x_obs_world);
				particle_item.sense_y.push_back(y_obs_world);
				particle_item.associations.push_back(1e6);   // Set association to land-mark to high-number initially
			}

		}  // End of looping over land-mark observation
		 
		// We have collect all qualified sensed detection for particular particular item
		// Now, we would like to associate which landmark that this detection correspond to
		AssociateObsLandmark(particle_item, map_landmarks, std_landmark);


		// Debug portion 
		// ==================================
	//	cout << "Particle ID : " << particle_item.id << "Weight = " << particle_item.weight << endl;
		// =================================

	}  // End of looping over particle


	// Normalize particle weight
	double weight_sum = 0.0;
	for (int i = 0; i < particles.size(); ++i)   // Get total sum of particle weigth
	{
		weight_sum = weight_sum + particles[i].weight;

	}
	
	// Debug only
	if (weight_sum == 0)
	{
		cout << " Error in data-associations " << endl;
		cout << "=============================" << endl;
		cout << "Particle 5 " << endl;

		cout << endl;
		cout << "Qualified landmark " << endl;

		int lmark_iter = 0;
		for (int a = 0; a < particles[5].qualified_lmark.size(); ++a)
		{
			lmark_iter = particles[5].qualified_lmark[a].first - 1;
			cout << "landmark id : " << map_landmarks.landmark_list[lmark_iter].id_i
				<< " , landmark xpos = " << map_landmarks.landmark_list[lmark_iter].x_f
				<< " , landmark ypos = " << map_landmarks.landmark_list[lmark_iter].y_f << endl;
		}

		cout << "Sense X" << "\t" << "Sense Y" << "\t" << "Lmark ID \n" ;
		for (int a = 0; a < particles[5].sense_x.size(); ++a)
		{
			cout << particles[5].sense_x[a] << "\t" << particles[5].sense_y[a] << "\t" << particles[5].associations[a] << endl;


		}



	}


	for (auto &particle_item : particles)          // Divide each particle weight with its total weight (normalize it)
	{
		particle_item.weight = particle_item.weight / weight_sum;

	}


}  // End of updateWeight function

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution

	int particle_index;   // Index generated for selecting particle during resampling process	
	struct temp_position_s {
		double x;
		double y;
		double theta;
	};                        // Temporary data-structure for saving position of particle, 
	
	temp_position_s dummy_var;

	// Define temporary data-structure for storing particle-filter position
	// during resampling process.
	// By the end of resampling, content of this temporary storage (x,y,theta) will be
	// copied back into original particles position. 
	std::vector<temp_position_s> temp_pos_storage;   // Temporary storage for particle position
	temp_pos_storage.reserve(particles.size());

	// Defining weight vector for resampling from particle-filter
	std::default_random_engine generator;
	std::vector<double> weight_vector;
	weight_vector.reserve(particles.size());


	// Populate weight vector, with weight from particle filter
	for (auto &particle_item : particles)
	{
		weight_vector.push_back(particle_item.weight);
	}
	// Initialize discrete-distribution machine with content of weight-vector
	std::discrete_distribution<> prob_weight(weight_vector.begin(), weight_vector.end());



	// Radomly sample from particle population, according to its weight 
	for (int i = 0; i < num_particles; ++i)
	{
		particle_index = prob_weight(generator);   // Index of particles

		temp_pos_storage.push_back(dummy_var);

		// Copy content of particles position specified by particle_index into temp_pos_storage
		temp_pos_storage[i].x = particles[particle_index].x;
		temp_pos_storage[i].y = particles[particle_index].y;
		temp_pos_storage[i].theta = particles[particle_index].theta;
	}


	// When re-sampling process is completed, copy back x,y,theta from temp_pos_storage into particles
	// position
	for (int i = 0; i < num_particles; ++i)
	{
		particles[i].x = temp_pos_storage[i].x;
		particles[i].y = temp_pos_storage[i].y;
		particles[i].theta = temp_pos_storage[i].theta;
	}


}

Particle ParticleFilter::SetAssociations(Particle& particle, const std::vector<int>& associations, 
                                     const std::vector<double>& sense_x, const std::vector<double>& sense_y)
{
    //particle: the particle to assign each listed association, and association's (x,y) world coordinates mapping to
    // associations: The landmark id that goes along with each listed association
    // sense_x: the associations x mapping already converted to world coordinates
    // sense_y: the associations y mapping already converted to world coordinates

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
