#include "sphere.h"
#include <stdlib.h>
#include <math.h>
#include <stdio.h>

extern Chessboard board;

/**********************************************************************
 * This function intersects a ray with a given sphere 'sph'. You should
 * use the parametric representation of a line and do the intersection.
 * The function should return the parameter value for the intersection,
 * which will be compared with others to determine which intersection
 * is closest. The value -1.0 is returned if there is no intersection
 *
 * If there is an intersection, the point of intersection should be
 * stored in the "hit" variable
 **********************************************************************/
float intersect_sphere(Point o, Vector u, Spheres *sph, Point *hit, bool inside) {
	// use the formula distance = (-b(+-)sqrt(b^2-4ac))/(2a)
	// a = u*u; b = 2*(l*(o-c)); c = (o-c)^2 - r^2

	Vector ctoo = get_vec(sph->center, o);
	float a = vec_dot(u, u);
	float b = 2 * vec_dot(u, ctoo);
	float c = vec_dot(ctoo, ctoo) - (sph->radius)*(sph->radius);

	float delta = b*b - 4*a*c;

	if(delta < 0) {
		// no intersection
		return -1.0;
	}	else {
		// has intersection, always pick the near one
		// if distance < 0(in the back), it will not be counted in intersect_scene()
		float distance;
		if(inside){
			// inside for refraction, take futher one
			distance = (-b + sqrt(delta)) / (2 * a);
		} else{
			distance = (-b - sqrt(delta)) / (2 * a);
		}
		*hit = get_point(o, vec_scale(u, distance));
		return distance;
	}
}

int board_colortest(Point hit){
	// return if the point hit is outside or white or black
	int x = (int)floor((hit.x - board.center.x)/board.grid_length);
	int z = (int)floor((hit.z - board.center.z)/board.grid_length);
	if (x >= 4 || x <= -5 || z >= 4 || z <= -5){
			return BOARD_OUTSIDE;
	} else if ((x+z) % 2 == 0) {
			return BOARD_BLACK;
	} else {
			return BOARD_WHITE;
	}
}

float intersect_board(Point o, Vector u, Point *hit) {
	// intersect with board
	Vector normal = board.surf_norm;
	if(vec_dot(normal, u) >= 0){
		// can't see board
		return -1.0;
	} else {
		float distance = (board.center.y-o.y) / u.y;
		Point hitBoard = get_point(o, vec_scale(u, distance));
		if (board_colortest(hitBoard) == BOARD_OUTSIDE){
				return -1.0;
		} else {
				*hit = hitBoard;
				return distance;
		}
	}
}

/*********************************************************************
 * This function returns a pointer to the sphere object that the
 * ray intersects first; NULL if no intersection. You should decide
 * which arguments to use for the function. For exmaple, note that you
 * should return the point of intersection to the calling function.
 **********************************************************************/
Spheres *intersect_scene(Point o, Vector u, Spheres *firstSphere, Point *hit, bool inside) {
	float min_d = 10000; // big enough
	Spheres* currentSphere = firstSphere;
	Spheres* resultSphere = NULL;

	while(currentSphere != NULL){
		Point tempHit;
		float current_d = -1.0;
		if(currentSphere->index == 0){
			current_d = intersect_board(o, u, &tempHit);
		} else{
			current_d = intersect_sphere(o, u, currentSphere, &tempHit, inside);
		}
		if(current_d > 0 && current_d < min_d) {
			// set thershold at 0 so the intersect point won't be point o
			min_d = current_d;
			*hit = tempHit;
			resultSphere = currentSphere;
		}
		currentSphere = currentSphere->next;
	}

	return resultSphere;
}

/*****************************************************
 * This function adds a sphere into the sphere list
 *
 * You need not change this.
 *****************************************************/
Spheres *add_sphere(Spheres *slist, Point ctr, float rad, float amb[],
		    float dif[], float spe[], float shine,
		    float refl, int sindex) {
  Spheres *new_sphere;

  new_sphere = (Spheres *)malloc(sizeof(Spheres));
  new_sphere->index = sindex;
  new_sphere->center = ctr;
  new_sphere->radius = rad;
  (new_sphere->mat_ambient)[0] = amb[0];
  (new_sphere->mat_ambient)[1] = amb[1];
  (new_sphere->mat_ambient)[2] = amb[2];
  (new_sphere->mat_diffuse)[0] = dif[0];
  (new_sphere->mat_diffuse)[1] = dif[1];
  (new_sphere->mat_diffuse)[2] = dif[2];
  (new_sphere->mat_specular)[0] = spe[0];
  (new_sphere->mat_specular)[1] = spe[1];
  (new_sphere->mat_specular)[2] = spe[2];
  new_sphere->mat_shineness = shine;
  new_sphere->reflectance = refl;
  new_sphere->next = NULL;

  if (slist == NULL) { // first object
    slist = new_sphere;
  } else { // insert at the beginning
    new_sphere->next = slist;
    slist = new_sphere;
  }

  return slist;
}

/******************************************
 * computes a sphere normal - done for you
 ******************************************/
Vector sphere_normal(Point q, Spheres *sph) {
  Vector rc;

  rc = get_vec(sph->center, q);
  normalize(&rc);
  return rc;
}

//--------------------------------------------------
// for adding a chess board, we fixes all parameters
Spheres *add_chessboard(Spheres *slist) {
		board.center = {0, -2.5, -4};
    board.surf_norm = {0, 1, 0};
    board.grid_length = 1.5;

		for(int i=0; i<3; i++){
			board.black_ambient[i] = 1.0;
		  board.black_diffuse[i] = 1.0;
		  board.black_specular[i] = 1.0;
		  board.white_ambient[i] = 0.0;
		  board.white_diffuse[i] = 0.0;
		  board.white_specular[i] = 0.0;
		}


    board.mat_shineness = 15;
    board.reflectance = 0.3;

    // add a fake sphere to sphere list with index 0 identifying chessboard
    slist = add_sphere(slist, board.center, 0.0, board.black_ambient, board.black_diffuse, board.black_specular, board.mat_shineness, board.reflectance, 0);
		return slist;
}
