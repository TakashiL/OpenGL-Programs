/**********************************************************************
 * Some stuff to handle spheres
 **********************************************************************/
#include "vector.h"

typedef struct sphere {
  int index;               // identifies a sphere; must be greater than 0

  Point center;
  float radius;

  float mat_ambient[3];    // material property used in Phong model
  float mat_diffuse[3];
  float mat_specular[3];
  float mat_shineness;

  float reflectance;
  float transparency = 1.5;
  float refractance;

  struct sphere *next;
} Spheres;   // a list of spheres

// intersect ray with sphere
Spheres *intersect_scene(Point, Vector, Spheres *, Point *, bool);
// return the unit normal at a point on sphere
Vector sphere_normal(Point, Spheres *);
// add a sphere to the sphere list
Spheres *add_sphere(Spheres *, Point, float, float [], float [], float [], float, float, int);

// --------------------------------------------

#define BOARD_BLACK 0
#define BOARD_WHITE 1
#define BOARD_OUTSIDE 2

typedef struct chessboard {
  Point center;
  Vector surf_norm;
  float grid_length;

  float black_ambient[3];
  float black_diffuse[3];
  float black_specular[3];
  float white_ambient[3];
  float white_diffuse[3];
  float white_specular[3];
  float mat_shineness;

  float reflectance;
} Chessboard; // one chessboard

int board_colortest(Point);
Spheres *add_chessboard(Spheres *);

// -----------------------------------------------

typedef struct face {
  Point p1;
  Point p2;
  Point p3;

} Faces; // a list of Faces (triangle mesh)

Spheres *add_face(Spheres *slist, Point p1, Point p2, Point p3, int index);
Spheres *read_chesspiece(Spheres *slist, Point, int pieceNum, const char *fileName);
Faces *findFace(int index);
Vector face_normal(Faces *face);
