#include "sphere.h"
#include "global.h"
#include <stdlib.h>
#include <math.h>
#include <stdio.h>

extern Chessboard board;
extern Faces faces[2][FACE_MAX_NUM];

Faces *findFace(int index);
Vector face_normal(Faces *face);

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
	if ((x+z) % 2 == 0) {
			return BOARD_WHITE;
	} else {
			return BOARD_BLACK;
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
		*hit = hitBoard;
		return distance;
	}
}

float intersect_face(Point o, Vector u, Faces *face, Point *hit){
	Vector normal = face_normal(face);
	if(vec_dot(normal, u) >= 0){
		// can't see face
		return -1.0;
	} else {
		Point p1 = face->p1;
		Point p2 = face->p2;
		Point p3 = face->p3;
		float distance = 1.0 * vec_dot(get_vec(o, p1), normal) / vec_dot(u, normal);
		Point hitface = get_point(o, vec_scale(u, distance));

		Vector w = get_vec(p1, hitface);
		Vector v = get_vec(p1, p2);
		Vector u = get_vec(p1, p3);

		float base = vec_dot(u,v)*vec_dot(u,v) - vec_dot(v,v)*vec_dot(u,u);
		float s = (vec_dot(u,v)*vec_dot(w,v) - vec_dot(v,v)*vec_dot(w,u)) / base;
		float t = (vec_dot(u,v)*vec_dot(w,u) - vec_dot(u,u)*vec_dot(w,v)) / base;

		if(s<0 || t<0 || s+t>1){
			return -1.0;
		}

		*hit = hitface;
		return distance;
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
			// intersect with face
			Faces *currentFace = findFace(currentSphere->index);
			current_d = intersect_face(o, u, currentFace, &tempHit);

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

Vector face_normal(Faces *face){
	Vector u = get_vec(face->p1, face->p2);
	Vector v = get_vec(face->p1, face->p3);

	float i = u.y*v.z - u.z*v.y;
	float j = u.z*v.x - u.x*v.z;
	float k = u.x*v.y - u.y*v.x;

	Vector rc = {i, j, k};
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

    board.mat_shineness = 40;
    board.reflectance = 0.5;

    // add a fake sphere to sphere list with index 0 identifying chessboard
    slist = add_sphere(slist, board.center, 0.0, board.black_ambient, board.black_diffuse, board.black_specular, board.mat_shineness, board.reflectance, 0);
		return slist;
}

//-------------------------------------------------
int pieceNum(int index){
	if(index/10000 == 1){
		// first piece
		return 0;
	} else{
		// second piece
		return 1;
	}
}

Faces *findFace(int index){
	return &faces[pieceNum(index)][index%10000];
}

Spheres *add_face(Spheres *slist, Point p1, Point p2, Point p3, int index){
	Faces *new_face = findFace(index);

	new_face->p1 = p1;
	new_face->p2 = p2;
	new_face->p3 = p3;

	float face_ambient[] = {0.7, 0.7, 0.7};
  float face_diffuse[] = {0.1, 0.5, 0.8};
  float face_specular[] = {1.0, 1.0, 1.0};
  float face_shineness = 10;
  float face_reflectance = 0.4;

  slist = add_sphere(slist, p1, 0.0, face_ambient, face_diffuse, face_specular, face_shineness, face_reflectance, index);
	slist->transparency = 1.5;
	slist->reflectance = 0.5;

	return slist;
}

Spheres *read_chesspiece(Spheres *slist, Point offset, int pieceNum, const char *fileName) {
	// pieceNum should be 0 or 1 identifying first or second piece
	FILE *dataFile;
	dataFile = fopen(fileName, "r");

	char identifier;
	int pointNum, faceNum;
	// read first line
	fscanf(dataFile, "%c %d %d ", &identifier, &pointNum, &faceNum);

	Point* points;
	points = new Point[pointNum];

	// read and store points information
	for(int i = 0; i < pointNum; i++) {
		float x, y, z;
		fscanf(dataFile, "%c %f %f %f ", &identifier, &x, &y, &z);
		points[i] = {x+offset.x, y+offset.y, z+offset.z};
	}

	// Read and store face informaton
	for(int i = 0; i < faceNum; i++) {
		int x, y, z;
		fscanf(dataFile, "%c %d %d %d ", &identifier, &x, &y, &z);
		int face_index = (pieceNum+1)*10000+i;
		slist = add_face(slist, points[x-1], points[y-1], points[z-1], face_index);
	}

	return slist;
}
