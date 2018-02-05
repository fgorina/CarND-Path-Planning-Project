#include <fstream>
#include <math.h>
#include <uWS/uWS.h>
#include <chrono>
#include <iostream>
#include <thread>
#include <vector>
#include "Eigen-3.3/Eigen/Core"
#include "Eigen-3.3/Eigen/QR"
#include "json.hpp"
#include "spline.h"

using namespace std;

// for convenience
using json = nlohmann::json;

// For converting back and forth between radians and degrees.
constexpr double pi() { return M_PI; }
double deg2rad(double x) { return x * pi() / 180; }
double rad2deg(double x) { return x * 180 / pi(); }

// Checks if the SocketIO event has JSON data.
// If there is data the JSON object in string format will be returned,
// else the empty string "" will be returned.
string hasData(string s) {
  auto found_null = s.find("null");
  auto b1 = s.find_first_of("[");
  auto b2 = s.find_first_of("}");
  if (found_null != string::npos) {
    return "";
  } else if (b1 != string::npos && b2 != string::npos) {
    return s.substr(b1, b2 - b1 + 2);
  }
  return "";
}

double distance(double x1, double y1, double x2, double y2)
{
	return sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1));
}
int ClosestWaypoint(double x, double y, const vector<double> &maps_x, const vector<double> &maps_y)
{

	double closestLen = 100000; //large number
	int closestWaypoint = 0;

	for(int i = 0; i < maps_x.size(); i++)
	{
		double map_x = maps_x[i];
		double map_y = maps_y[i];
		double dist = distance(x,y,map_x,map_y);
		if(dist < closestLen)
		{
			closestLen = dist;
			closestWaypoint = i;
		}

	}

	return closestWaypoint;

}

int NextWaypoint(double x, double y, double theta, const vector<double> &maps_x, const vector<double> &maps_y)
{

	int closestWaypoint = ClosestWaypoint(x,y,maps_x,maps_y);

	double map_x = maps_x[closestWaypoint];
	double map_y = maps_y[closestWaypoint];

	double heading = atan2((map_y-y),(map_x-x));

	double angle = fabs(theta-heading);
  angle = min(2*pi() - angle, angle);

  if(angle > pi()/4)
  {
    closestWaypoint++;
  if (closestWaypoint == maps_x.size())
  {
    closestWaypoint = 0;
  }
  }

  return closestWaypoint;
}

// Transform from Cartesian x,y coordinates to Frenet s,d coordinates
vector<double> getFrenet(double x, double y, double theta, const vector<double> &maps_x, const vector<double> &maps_y)
{
	int next_wp = NextWaypoint(x,y, theta, maps_x,maps_y);

	int prev_wp;
	prev_wp = next_wp-1;
	if(next_wp == 0)
	{
		prev_wp  = maps_x.size()-1;
	}

	double n_x = maps_x[next_wp]-maps_x[prev_wp];
	double n_y = maps_y[next_wp]-maps_y[prev_wp];
	double x_x = x - maps_x[prev_wp];
	double x_y = y - maps_y[prev_wp];

	// find the projection of x onto n
	double proj_norm = (x_x*n_x+x_y*n_y)/(n_x*n_x+n_y*n_y);
	double proj_x = proj_norm*n_x;
	double proj_y = proj_norm*n_y;

	double frenet_d = distance(x_x,x_y,proj_x,proj_y);

	//see if d value is positive or negative by comparing it to a center point

	double center_x = 1000-maps_x[prev_wp];
	double center_y = 2000-maps_y[prev_wp];
	double centerToPos = distance(center_x,center_y,x_x,x_y);
	double centerToRef = distance(center_x,center_y,proj_x,proj_y);

	if(centerToPos <= centerToRef)
	{
		frenet_d *= -1;
	}

	// calculate s value
	double frenet_s = 0;
	for(int i = 0; i < prev_wp; i++)
	{
		frenet_s += distance(maps_x[i],maps_y[i],maps_x[i+1],maps_y[i+1]);
	}

	frenet_s += distance(0,0,proj_x,proj_y);

	return {frenet_s,frenet_d};

}

// Transform from Frenet s,d coordinates to Cartesian x,y
vector<double> getXY(double s, double d, const vector<double> &maps_s, const vector<double> &maps_x, const vector<double> &maps_y)
{
	int prev_wp = -1;

	while(s > maps_s[prev_wp+1] && (prev_wp < (int)(maps_s.size()-1) ))
	{
		prev_wp++;
	}

	int wp2 = (prev_wp+1)%maps_x.size();

	double heading = atan2((maps_y[wp2]-maps_y[prev_wp]),(maps_x[wp2]-maps_x[prev_wp]));
	// the x,y,s along the segment
	double seg_s = (s-maps_s[prev_wp]);

	double seg_x = maps_x[prev_wp]+seg_s*cos(heading);
	double seg_y = maps_y[prev_wp]+seg_s*sin(heading);

	double perp_heading = heading-pi()/2;

	double x = seg_x + d*cos(perp_heading);
	double y = seg_y + d*sin(perp_heading);

	return {x,y};

}

vector<vector<double>> generate_trajectory(double car_x, double car_y, double car_s, double car_d, double car_yaw,
    vector<double> previous_path_x, vector<double>previous_path_y, double pts_space,
    int from_lane, int to_lane, double ref_vel, double target_x, int n_points, double delta,
     const vector<double> &map_waypoints_s, const vector<double> &map_waypoints_x, const vector<double> &map_waypoints_y){

    int prev_size = previous_path_x.size();
    // 30m spaced points in the desired pathway

    vector<double> ptsx;
    vector<double> ptsy;

    // Initial car position

    double ref_x = car_x;
    double ref_y = car_y;
    double ref_yaw = deg2rad(car_yaw);

    if (prev_size < 2){
        double prev_car_x = ref_x - cos(car_yaw);
        double prev_car_y = ref_y - sin(car_yaw);

        ptsx.push_back(prev_car_x);
        ptsx.push_back(car_x);
        ptsy.push_back(prev_car_y);
        ptsy.push_back(car_y);
    }else {
        ref_x = previous_path_x[prev_size-1];
        ref_y = previous_path_y[prev_size-1];

        double ref_x_prev = previous_path_x[prev_size-2];
        double ref_y_prev = previous_path_y[prev_size-2];
        ref_yaw = atan2(ref_y-ref_y_prev, ref_x - ref_x_prev);


        ptsx.push_back(ref_x_prev);
        ptsx.push_back(ref_x);
        ptsy.push_back(ref_y_prev);
        ptsy.push_back(ref_y);

    }


    // Generate trajectory

    vector<double> next_x_vals;
    vector<double> next_y_vals;

    for(int i = 0; i < 3; i++){

        vector<double> next_point = getXY(car_s+pts_space*(i+1), (2+4*to_lane), map_waypoints_s, map_waypoints_x, map_waypoints_y);

        ptsx.push_back(next_point[0]);
        ptsy.push_back(next_point[1]);
    }

    for(int i = 0; i < ptsx.size(); i++){

        double shift_x = ptsx[i] - ref_x;
        double shift_y = ptsy[i] - ref_y;

        ptsx[i] = (shift_x * cos(0 - ref_yaw) - shift_y * sin(0 - ref_yaw));
        ptsy[i] = (shift_x * sin(0 - ref_yaw) + shift_y * cos(0 - ref_yaw));
    }

    tk::spline s;

    /**
    cout << "--------------------------" << endl;
    for(int i = 0; i < ptsx.size(); i++){
        cout << "x " << ptsx[i] << endl;
    }
    */

    s.set_points(ptsx, ptsy);

    for(int i = 0; i < previous_path_x.size(); i++){
        next_x_vals.push_back(previous_path_x[i]);
        next_y_vals.push_back(previous_path_y[i]);
    }

    double target_y = s(target_x);
    double target_dist = sqrt((target_x)*(target_x)+(target_y)*(target_y));

    double x_add_on = 0;


    double N = (target_dist)/(delta*ref_vel);

    //        cout << "N " << N << " delta v" << delta_v << " speed " << car_speed << " new speed " << ref_vel << endl;

    for(int i = 1; i <= n_points-previous_path_x.size(); i++){
        double x_point = x_add_on + (target_x) / N;
        double y_point = s(x_point);

        x_add_on = x_point;

        double x_ref = x_point;
        double y_ref = y_point;

        x_point = (x_ref * cos(ref_yaw) - y_ref * sin(ref_yaw));
        y_point = (x_ref * sin(ref_yaw) + y_ref * cos(ref_yaw));

        x_point += ref_x;
        y_point += ref_y;
        next_x_vals.push_back(x_point);
        next_y_vals.push_back(y_point);

    }

    vector<vector<double>> out;
    out.push_back(next_x_vals);
    out.push_back(next_y_vals);
    return(out);
}

int main() {
  uWS::Hub h;

  // Load up map values for waypoint's x,y,s and d normalized normal vectors
  vector<double> map_waypoints_x;
  vector<double> map_waypoints_y;
  vector<double> map_waypoints_s;
  vector<double> map_waypoints_dx;
  vector<double> map_waypoints_dy;

  // Waypoint map to read from
  string map_file_ = "../data/highway_map.csv";
  // The max s value before wrapping around the track back to 0
  double max_s = 6945.554;

  ifstream in_map_(map_file_.c_str(), ifstream::in);

  string line;
  while (getline(in_map_, line)) {
  	istringstream iss(line);
  	double x;
  	double y;
  	float s;
  	float d_x;
  	float d_y;
  	iss >> x;
  	iss >> y;
  	iss >> s;
  	iss >> d_x;
  	iss >> d_y;
  	map_waypoints_x.push_back(x);
  	map_waypoints_y.push_back(y);
  	map_waypoints_s.push_back(s);
  	map_waypoints_dx.push_back(d_x);
  	map_waypoints_dy.push_back(d_y);
  }


  int lane = 1; // Desired lane
  double ref_vel = 0.0;    // Desired speed in this  timeframe
  int timestep = 0;

  double acum_vd = 0.0;
  double acum_vd2 = 0.0;
  int npoints = 0.0;
  int state = 0; // 0-> free driving, 1-> following car 2-> slowing for distance, 3->changinh right


  h.onMessage([&map_waypoints_x,&map_waypoints_y,&map_waypoints_s,&map_waypoints_dx,&map_waypoints_dy,
                &lane, &ref_vel, &timestep, &acum_vd, &acum_vd2, &npoints, &state]
                (uWS::WebSocket<uWS::SERVER> ws, char *data, size_t length,
                     uWS::OpCode opCode) {
    // "42" at the start of the message means there's a websocket message event.
    // The 4 signifies a websocket message
    // The 2 signifies a websocket event
    //auto sdata = string(data).substr(0, length);
    //cout << sdata << endl;
    if (length && length > 2 && data[0] == '4' && data[1] == '2') {

      auto s = hasData(data);

      if (s != "") {
        auto j = json::parse(s);

        string event = j[0].get<string>();

        if (event == "telemetry") {
          // j[1] is the data JSON object

        	// Main car's localization Data
          	double car_x = j[1]["x"];
          	double car_y = j[1]["y"];
          	double car_s = j[1]["s"];
          	double car_d = j[1]["d"];
          	double car_yaw = j[1]["yaw"];
          	double car_speed = j[1]["speed"];
            car_speed = car_speed / 2.24; // Everything in m/s

          	// Previous path data given to the Planner
          	auto previous_path_x_orig = j[1]["previous_path_x"];
          	auto previous_path_y_orig = j[1]["previous_path_y"];
          	// Previous path's end s and d values
          	double end_path_s = j[1]["end_path_s"];
          	double end_path_d = j[1]["end_path_d"];

          	// Sensor Fusion Data, a list of all other cars on the same side of the road.
          	auto sensor_fusion = j[1]["sensor_fusion"];

          	json msgJson;

            vector<double> previous_path_x;
            vector<double> previous_path_y;

            int prev_path_max_size = 5;
            double max_vel = 49.5/ 2.24 ;  // Desired maximum speed
            double delta = 0.02;        // simulator delta t
            double max_forward_acc = 5.0 ; // en m/s
            double max_brake_acc = 8.0; // en m/s sembla ok
            double pts_space = 30.0;  // punts de referencia espaiats 30 m
            double next_car_distance = 15.0;    // Minima distància al cotxe del dvant .. 15 funciona 20 llarg
            double prev_car_distance = 15.0;    // Mínima distànciadel cotxe previ. Eren 7
            double speed_horizon = 30.0;    // mtrs. at which to look for fastest lane. Sensors?
            int next_lane = lane;       // By default follow same lane
            bool europe = true; // Then try to drive to the rightest lane and NEVER advance by the right

            // Copy first 5 items if they next_car_distance

            for(int i = 0; i < previous_path_x_orig.size() && i < prev_path_max_size; i++){
                double xi = previous_path_x_orig[i];
                double yi = previous_path_y_orig[i];
                previous_path_x.push_back(xi);
                previous_path_y.push_back(yi);
            }
          	// TODO: define a path made up of (x,y) points that the car will visit sequentially every .02 seconds

            timestep += 1;  // Increase timestep

            int prev_size = previous_path_x.size(); // Alrady passed path point

            // Treat Sensor sensor_fusion


            if (prev_size > 0){
                car_s = end_path_s;
            }

            bool too_close = false;

            vector<bool> free_lane;         // true if there is place to change lanes. Initilized to true
            free_lane.push_back(true);
            free_lane.push_back(true);
            free_lane.push_back(true);

            vector<double> speed_lane;      // speed of next car. Initialized to max_vel
            speed_lane.push_back(max_vel);
            speed_lane.push_back(max_vel);
            speed_lane.push_back(max_vel);

            // Check free lanes and compute lane speed

            for (int i = 0; i < sensor_fusion.size(); i++){

                float d = sensor_fusion[i][6];  // d coordinate of other car

                // Compute lane of car

                int car_lane = floor(d / 4);

                if (car_lane >= 0 && car_lane < 3){

                    double vx = sensor_fusion[i][3];
                    double vy = sensor_fusion[i][4];
                    double cs = sensor_fusion[i][5];
                    double cd = sensor_fusion[i][6];

                    double check_speed = sqrt(vx*vx+vy*vy);

                    double check_car_s = sensor_fusion[i][5];
                    check_car_s += ((double)prev_size * delta * check_speed);


                    /** Uncomment to get stats about lateral speeds
                    // Compute vs and vd. First

                    vector<double> orig = getXY(cs, cd, map_waypoints_s, map_waypoints_x, map_waypoints_y);
                    vector<double> nexts = getXY(cs+1.0, cd, map_waypoints_s, map_waypoints_x, map_waypoints_y);
                    vector<double> nextd = getXY(cs, cd+1.0, map_waypoints_s, map_waypoints_x, map_waypoints_y);

                    // Unit vectors in s and d direction

                    double usx = nexts[0]-orig[0];
                    double usy = nexts[1]-orig[1];
                    double udx = nextd[0]-orig[0];
                    double udy = nextd[1]-orig[1];

                    // Now we may project speed along d unit vector so we get vd

                    double vd = vx*udx + vy*udy;
                    npoints += 1;
                    acum_vd += vd;
                    acum_vd2 += vd*vd;
                    */

                    // A car lowers speed to all lanes which are "near" it so
                    // a car changing lanes lowers speed of 2 lanes, origin and destination
                    // It is detected because it is nearer thah 3m to a center lane
                    if((check_car_s > car_s) &&
                        ((check_car_s - car_s ) < (speed_horizon))){   // Far away cars doesnt mean nothing
                            for(int i = 0; i < 3; i++){
                                if (fabs(d - (2 + 4*i)) < 3){
                                    speed_lane[i] = min(speed_lane[i], check_speed);
                                }
                            }
                    }

                    double pr = prev_car_distance;

                    if (check_speed > car_speed){  // Car is faster than we are
                        pr = pr + (check_speed - car_speed) * 2.0; // 1.5 is a factor that comes from the time to run the planning distance (30m)
                    }
                    if((check_car_s < car_s ) && ((car_s - check_car_s) < pr)){
                        for(int i = 0; i < 3; i++){
                            if (fabs(d - (2 + 4*i)) < 3){
                                free_lane[i] = false;
                            }
                        }
                    }

                    if((check_car_s >= car_s ) && ((check_car_s - car_s) < next_car_distance)){
                        for(int i = 0; i < 3; i++){
                            if (fabs(d - (2 + 4*i)) < 3){
                                free_lane[i] = false;
                            }
                        }
                    }
                }
            }

            double target_speed = max_vel;
            double check_distance = next_car_distance;

            // If we are already following a car we need to be twice as Far
            // to release the state. It is like a histeresis.
            if(state == 1){
                check_distance = 2 * check_distance;
                state = 0;
            }

            for (int i = 0; i < sensor_fusion.size(); i++){

                float d = sensor_fusion[i][6];  // d coordinate of other car

                //if (d < (2 + (4*lane) + 2) && d > (2 + (4*lane) - 2)){     // My Lane

                if (fabs(car_d - d) <= 3.0){     // My Lane

                    double vx = sensor_fusion[i][3];
                    double vy = sensor_fusion[i][4];
                    double check_speed = sqrt(vx*vx+vy*vy) ; // Seems is m/s convert to m/h
                    double check_car_s = sensor_fusion[i][5];

                    check_car_s += ((double)prev_size * delta * check_speed);

                    // Here we detect if we have a car in front which goes slower or is too near
                    if ((check_car_s >= car_s) && (
                        (fabs(check_car_s - car_s) < check_distance && (check_speed <= car_speed)) || fabs(check_car_s - car_s) < next_car_distance)){
                        if ( fabs(check_car_s - car_s) < next_car_distance*0.8){
                            too_close = true;    // I would prefer put the other car speed
                        }
                        state = 1;
                        target_speed = min(target_speed, check_speed);

                        //cout << "Car at " << check_car_s - car_s << " m" << " his speed " << target_speed  << " my speed " << car_speed << endl;
                    }
                }
            }


            // Now we print average acum_vd and stdev acum_vd2

            //double avg_vd = acum_vd / double(npoints);
            //double sigma2_vd = (acum_vd2 + npoints*avg_vd*avg_vd-2.0*avg_vd*acum_vd)/npoints;

            //cout << "Average " << avg_vd << " sigma2 " << sigma2_vd << endl;

            // Here we compute which is the fastest lane

            double fastest_speed = speed_lane[lane];
            double fastest_lane = lane;

            // Try to get a faster lane
            for (int i = 0; i < 3; i++){
                if (speed_lane[i] > fastest_speed){
                    fastest_speed = speed_lane[i];
                    fastest_lane = i;
                }
            }

            //cout << "Fastest Lane " << fastest_lane << endl;

            if (state == 1 ){    // Reduce speed before changing lane.

                if (ref_vel < target_speed  && !too_close){
                    ref_vel += max_forward_acc * delta / 2.0;
                }else {
                    //cout << "Braking" << endl;
                    ref_vel -=  max_brake_acc * delta;
                }
            }

            // If in normal state just increase till max_vel

            if (state == 0 && ref_vel < max_vel){
                ref_vel += max_forward_acc * delta;
            }

            // Right comparison is because if not it may generate some too strong accelerations
            // when jumping 2 lanes, so we made 2 changes instead of just one

            if (fastest_lane > lane && free_lane[lane+1] && (car_d - (2 + lane*4)) > -0.5){
                next_lane = lane + 1;
                state = 0;
            }
            /** else if (fastest_lane > lane && !free_lane[lane+1] && state == 1){
                if (ref_vel > speed_lane[lane+1]*0.8){
                    ref_vel -=  max_brake_acc * delta;
                    cout << "Special braking " << endl;
                }
            }
            */
          else if (fastest_lane < lane && free_lane[lane-1] && ((2 + lane*4) - car_d ) > -0.5){
                next_lane = lane-1;
                state = 0;
            }

            /** else if (fastest_lane < lane && !free_lane[lane-1] && state == 1){
                if (ref_vel > speed_lane[lane-1]*0.8){
                    ref_vel -=  max_brake_acc * delta;
                    cout << "Special braking " << endl;
                }
            }
            */


            //  Block lane change for small difference speed

            if (fastest_speed < (speed_lane[lane]+0.1)){
               next_lane = lane;
            }


            if (lane != next_lane){
                cout <<  " t = " << timestep << " Changing from lane " << lane << " to lane " << next_lane << endl;
            }

            lane =  next_lane;

            // Generate trajectory

            vector<vector<double>> trj= generate_trajectory(car_x, car_y, car_s, car_d, car_yaw,
                previous_path_x, previous_path_y, pts_space,
                lane, next_lane, ref_vel, 20.0, 30, delta,
                 map_waypoints_s, map_waypoints_x, map_waypoints_y);



            msgJson["next_x"] = trj[0];
            msgJson["next_y"] = trj[1];

          	auto msg = "42[\"control\","+ msgJson.dump()+"]";

          	//this_thread::sleep_for(chrono::milliseconds(1000));
          	ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);

        }
      } else {
        // Manual driving
        std::string msg = "42[\"manual\",{}]";
        ws.send(msg.data(), msg.length(), uWS::OpCode::TEXT);
      }
    }
  });

  // We don't need this since we're not using HTTP but if it's removed the
  // program
  // doesn't compile :-(
  h.onHttpRequest([](uWS::HttpResponse *res, uWS::HttpRequest req, char *data,
                     size_t, size_t) {
    const std::string s = "<h1>Hello world!</h1>";
    if (req.getUrl().valueLength == 1) {
      res->end(s.data(), s.length());
    } else {
      // i guess this should be done more gracefully?
      res->end(nullptr, 0);
    }
  });

  h.onConnection([&h](uWS::WebSocket<uWS::SERVER> ws, uWS::HttpRequest req) {
    std::cout << "Connected!!!" << std::endl;
  });

  h.onDisconnection([&h](uWS::WebSocket<uWS::SERVER> ws, int code,
                         char *message, size_t length) {
    ws.close();
    std::cout << "Disconnected" << std::endl;
  });

  int port = 4567;
  if (h.listen(port)) {
    std::cout << "Listening to port " << port << std::endl;
  } else {
    std::cerr << "Failed to listen to port" << std::endl;
    return -1;
  }
  h.run();
}
