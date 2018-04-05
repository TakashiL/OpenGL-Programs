#include <stdio.h>
#include <GL/glut.h>
#include <math.h>
#include "global.h"
#include "sphere.h"

//
// Global variables
//
extern int win_width;
extern int win_height;

extern GLfloat frame[WIN_HEIGHT][WIN_WIDTH][3];

extern float image_width;
extern float image_height;

extern Point eye_pos;
extern float image_plane;
extern RGB_float background_clr;
extern RGB_float null_clr;

extern Spheres *scene;
extern Chessboard board;

// light 1 position and color
extern Point light1;
extern float light1_ambient[3];
extern float light1_diffuse[3];
extern float light1_specular[3];

// global ambient term
extern float global_ambient[3];

// light decay parameters
extern float decay_a;
extern float decay_b;
extern float decay_c;

extern int step_max;
extern int shadow_on;
extern int reflection_on;
extern int refraction_on;
extern int chess_board_on;
extern int stochastic_on;
extern int supersampling_on;

/////////////////////////////////////////////////////////////////////

/*********************************************************************
 * Phong illumination - q is the point on the sphere that we want to
 get the color, viewtoq is the vector from view point to q, surf_norm is n.*
 *********************************************************************/
RGB_float phong(Point q, Vector viewtoq, Vector surf_norm, Spheres *sph, bool inside) {
	RGB_float color;

	if (sph->index == 0) {
		// update matrices for chess board
		int board_color = board_colortest(q);
		for (int i = 0; i < 3; i++) {
			if(board_color == BOARD_BLACK){
				sph->mat_ambient[i] = board.black_ambient[i];
				sph->mat_diffuse[i] = board.black_diffuse[i];
				sph->mat_specular[i] = board.black_specular[i];
			} else{
				sph->mat_ambient[i] = board.white_ambient[i];
				sph->mat_diffuse[i] = board.white_diffuse[i];
				sph->mat_specular[i] = board.white_specular[i];
			}
		}
	}

	// normalize v and n anyway
	Vector v = vec_scale(viewtoq, -1.0); // q to view point
	normalize(&v);
	normalize(&surf_norm);

	// get l: vector from q to light source
	Vector l = get_vec(q, light1);
	normalize(&l);

	// to get r, use formula r = 2*(n*l)*n - l
	// r: vector of perfect reflection
	Vector r = vec_minus(vec_scale(surf_norm, 2.0*vec_dot(surf_norm, l)), l);
	normalize(&r);

	float d = vec_len(get_vec(q, light1));
	float decay_factor = 1.0/(decay_a + decay_b*d + decay_c*d*d);

	RGB_float Iga = {global_ambient[0], global_ambient[1], global_ambient[2]};
	RGB_float color_ga = clr_sepscale(Iga, sph->mat_ambient);

	RGB_float Ila = {light1_ambient[0], light1_ambient[1], light1_ambient[2]};
	RGB_float color_la = clr_sepscale(Ila, sph->mat_ambient);

	RGB_float Ild = {light1_diffuse[0], light1_diffuse[1], light1_diffuse[2]};
	RGB_float color_diffuse;
	// diffuse should be 0 if (n*l)<0
	if(vec_dot(surf_norm, l) < 0){
		color_diffuse = {0, 0, 0};
	} else {
		color_diffuse = clr_scale(clr_sepscale(Ild, sph->mat_diffuse), vec_dot(surf_norm, l));
	}

	RGB_float Ls = {light1_specular[0], light1_specular[1], light1_specular[2]};
	RGB_float color_specular;
	// specular should be 0 if (n*l)<0 or (r*v)<0
	if(vec_dot(surf_norm, l) < 0 || vec_dot(r, v) < 0){
		color_specular = {0, 0, 0};
	} else {
		color_specular = clr_scale(clr_sepscale(Ls, sph->mat_specular), pow(vec_dot(r, v), sph->mat_shineness));
	}

	RGB_float color_decaied = clr_scale(clr_add(color_diffuse, color_specular), decay_factor);

	// calculate final color
	color = clr_add(color_ga, color_la);

	// shadows
	bool in_shadow = false;
	if(shadow_on){
		Point blockPoint;
		Spheres *blockSphere = intersect_scene(q, l, scene, &blockPoint, inside); // this function will not return point q
		if(blockSphere != NULL){
			// in shadow
			in_shadow = true;
		}
	}
	if(!in_shadow && !inside){
		color = clr_add(color, color_decaied);
	}

	return color;
}

/************************************************************************
 * This is the recursive ray tracer - you need to implement this!
 * You should decide what arguments to use.
 ************************************************************************/
RGB_float recursive_ray_trace(Point eye, Vector ray, int step, bool inside) {
	RGB_float color;

	// find the intersect point
	Point hitPoint;
	Spheres *intersectSphere = intersect_scene(eye, ray, scene, &hitPoint, inside);

	if(intersectSphere == NULL){
		color = background_clr;
		return color;
	}

	// has intersection, get the color of intersect point
	Vector surf_norm;
	if (intersectSphere->index == 0) {
			//intersect with board
			surf_norm = board.surf_norm;
	} else {
			surf_norm = sphere_normal(hitPoint, intersectSphere);
			if(inside){
				surf_norm = vec_scale(surf_norm, -1.0);
			}
	}
	color = phong(hitPoint, ray, surf_norm, intersectSphere, inside);

	if(step >= step_max){
		return color;
	}

	// reflections
	if(reflection_on){
		Vector l = get_vec(hitPoint, eye);
		normalize(&l);
		Vector reflection_ray = vec_minus(vec_scale(surf_norm, 2.0*vec_dot(surf_norm, l)), l);
		normalize(&reflection_ray);
		RGB_float reflection_color = recursive_ray_trace(hitPoint, reflection_ray, step+1, inside);
		color = clr_add(color, clr_scale(reflection_color, intersectSphere->reflectance));
	}

	// refractions
	if(refraction_on && intersectSphere->index != 0){
		float snellRatio = 1.0 / intersectSphere->transparency;
		if(inside){
			snellRatio = 1.0 / snellRatio;
		}
		Vector refraction_i = ray;
		Vector refraction_n = surf_norm;
		normalize(&refraction_i);
		normalize(&refraction_n);

		float cos_theta = -1.0*vec_dot(refraction_i, refraction_n);
		float sin_theta_2 = snellRatio*snellRatio*(1-cos_theta*cos_theta);
		float para_n = snellRatio*cos_theta - sqrt(1-sin_theta_2);

		Vector refraction_ray = vec_plus(vec_scale(refraction_i, snellRatio), vec_scale(refraction_n, para_n));
		normalize(&refraction_ray);
		RGB_float refraction_color = recursive_ray_trace(hitPoint, refraction_ray, step+1, !(inside)); // change value of inside
		color = clr_add(color, clr_scale(refraction_color, intersectSphere->refractance));
	}

	// stochastic ray
	if(stochastic_on){
		RGB_float stoch_color_final = {0,0,0};
		int count = 0;
		while(count < STOCH_RAY_NUM){
			Vector stoch_ray = {rand()%10-5, rand()%10-5, rand()%10-5};
			while(vec_dot(stoch_ray, surf_norm) <= 0){
				stoch_ray = {rand()%10-5, rand()%10-5, rand()%10-5}; // renew the ray
			}
			normalize(&stoch_ray);
			RGB_float stoch_color = recursive_ray_trace(hitPoint, stoch_ray, step+1, inside); // change value of inside
			stoch_color_final = clr_add(stoch_color_final, clr_scale(stoch_color, intersectSphere->reflectance));
			count++;
		}
		color = clr_add(color, clr_scale(stoch_color_final, 1.0/STOCH_RAY_NUM));
	}

	return color;
}

/*********************************************************************
 * This function traverses all the pixels and cast rays. It calls the
 * recursive ray tracer and assign return color to frame
 *
 * You should not need to change it except for the call to the recursive
 * ray tracer. Feel free to change other parts of the function however,
 * if you must.
 *********************************************************************/
void ray_trace() {
  int i, j;
  float x_grid_size = image_width / float(win_width);
  float y_grid_size = image_height / float(win_height);
  float x_start = -0.5 * image_width;
  float y_start = -0.5 * image_height;
  RGB_float ret_color;
  Point cur_pixel_pos;
  Vector ray;

  // ray is cast through center of pixel
  cur_pixel_pos.x = x_start + 0.5 * x_grid_size;
  cur_pixel_pos.y = y_start + 0.5 * y_grid_size;
  cur_pixel_pos.z = image_plane;

  for (i=0; i<win_height; i++) {
    for (j=0; j<win_width; j++) {
      ray = get_vec(eye_pos, cur_pixel_pos);
      ret_color = recursive_ray_trace(eye_pos, ray, 0, false);

			if(supersampling_on){
				// add four more rays
				float delta_x = 0.25 * x_grid_size;
				float delta_y = 0.25 * y_grid_size;

				Point lu = {cur_pixel_pos.x - delta_x, cur_pixel_pos.y + delta_y, cur_pixel_pos.z};
				Vector lu_ray = get_vec(eye_pos, lu);
				ret_color = clr_add(ret_color, recursive_ray_trace(eye_pos, lu_ray, 0, false));

				Point lb = {cur_pixel_pos.x - delta_x, cur_pixel_pos.y - delta_y, cur_pixel_pos.z};
				Vector lb_ray = get_vec(eye_pos, lb);
				ret_color = clr_add(ret_color, recursive_ray_trace(eye_pos, lb_ray, 0, false));

				Point ru = {cur_pixel_pos.x + delta_x, cur_pixel_pos.y + delta_y, cur_pixel_pos.z};
				Vector ru_ray = get_vec(eye_pos, ru);
				ret_color = clr_add(ret_color, recursive_ray_trace(eye_pos, ru_ray, 0, false));

				Point rb = {cur_pixel_pos.x + delta_x, cur_pixel_pos.y - delta_y, cur_pixel_pos.z};
				Vector rb_ray = get_vec(eye_pos, rb);
				ret_color = clr_add(ret_color, recursive_ray_trace(eye_pos, rb_ray, 0, false));
			}


      // ret_color = background_clr; // just background for now

      // Parallel rays can be cast instead using below
      //
      // ray.x = ray.y = 0;
      // ray.z = -1.0;
      // ret_color = recursive_ray_trace(cur_pixel_pos, ray, 1);

			// Checkboard for testing
			// RGB_float clr = {float(i/32), 0, float(j/32)};
			// ret_color = clr;

      frame[i][j][0] = GLfloat(ret_color.r);
      frame[i][j][1] = GLfloat(ret_color.g);
      frame[i][j][2] = GLfloat(ret_color.b);

      cur_pixel_pos.x += x_grid_size;
    }

    cur_pixel_pos.y += y_grid_size;
    cur_pixel_pos.x = x_start;
  }
}
