
#include "Angel.h"
#include <string>

typedef Angel::vec4 point4;
typedef Angel::vec4 color4;

const int NumRobot = 36; //(6 faces)(2 triangles/face)(3 vertices/triangle)
const int NumSphereLoop = 342;
const int NumSphereFan = 40;

point4 robot_points[NumRobot];
color4 robot_colors[NumRobot];

point4 sphereLoop_points[NumSphereLoop];
color4 sphereLoop_colors[NumSphereLoop];

point4 sphereFan_points[NumSphereFan];
color4 sphereFan_colors[NumSphereFan];

point4 vertices[8] = {
    point4( -0.5, -0.5, 0.5, 1.0),
    point4( -0.5,  0.5, 0.5, 1.0),
    point4(  0.5,  0.5, 0.5, 1.0),
    point4(  0.5, -0.5, 0.5, 1.0),
    point4( -0.5, -0.5, -0.5, 1.0),
    point4( -0.5,  0.5, -0.5, 1.0),
    point4(  0.5,  0.5, -0.5, 1.0),
    point4(  0.5, -0.5, -0.5, 1.0)
};

// RGBA olors
color4 vertex_colors[8] = {
    color4( 0.0, 0.0, 0.0, 1.0),  // black
    color4( 1.0, 0.0, 0.0, 1.0),  // red
    color4( 1.0, 1.0, 0.0, 1.0),  // yellow
    color4( 0.0, 1.0, 0.0, 1.0),  // green
    color4( 0.0, 0.0, 1.0, 1.0),  // blue
    color4( 1.0, 0.0, 1.0, 1.0),  // magenta
    color4( 1.0, 0.6, 0.6, 1.0),  // pink
    color4( 0.0, 1.0, 1.0, 1.0)   // cyan
};


// Parameters controlling the size of the Robot's arm
const GLfloat BASE_HEIGHT      = 2.0;
const GLfloat BASE_WIDTH       = 5.0;
const GLfloat LOWER_ARM_HEIGHT = 5.0;
const GLfloat LOWER_ARM_WIDTH  = 0.5;
const GLfloat UPPER_ARM_HEIGHT = 5.0;
const GLfloat UPPER_ARM_WIDTH  = 0.5;

// Shader transformation matrices
mat4 model_view;
mat4 sphere_view;
mat4 final_sphereview;
GLuint ModelView, Projection;

// Array of rotation angles (in degrees) for each rotation axis
enum { Base = 0, LowerArm = 1, UpperArm = 2, NumAngles = 3 };
int Axis = Base;
GLfloat Theta[NumAngles] = {0.0};

// global Parameters for sphere
float old_x;
float old_y;
float old_z;
float new_x;
float new_y;
float new_z;
GLfloat fetchBaseAngle;
GLfloat fetchLowerAngle;
GLfloat fetchUpperAngle;
GLfloat putBaseAngle;
GLfloat putLowerAngle;
GLfloat putUpperAngle;
int basefetch_count = 0;
int lowerfetch_count = 0;
int upperfetch_count = 0;
int baseput_count = 0;
int lowerput_count = 0;
int upperput_count = 0;

bool topview_flag = false; // false means side view, true means top view
int status; // 1:fetch, 2:put, 3:go back
float speed = 10.0;

// location of vertex attributes in the shader program
GLuint vPosition;
GLuint vColor;

// VectorArrayObject and VectorBufferObject
GLuint vaos[3]; // Robot for the first and Sphere for the later two
GLuint vbos[3]; // Position and color for each vao

//----------------------------------------------------------------------------

int Index = 0;

void quad(int a, int b, int c, int d) {
    robot_colors[Index] = vertex_colors[a]; robot_points[Index] = vertices[a]; Index++;
    robot_colors[Index] = vertex_colors[a]; robot_points[Index] = vertices[b]; Index++;
    robot_colors[Index] = vertex_colors[a]; robot_points[Index] = vertices[c]; Index++;
    robot_colors[Index] = vertex_colors[a]; robot_points[Index] = vertices[a]; Index++;
    robot_colors[Index] = vertex_colors[a]; robot_points[Index] = vertices[c]; Index++;
    robot_colors[Index] = vertex_colors[a]; robot_points[Index] = vertices[d]; Index++;
}

void colorcube() {
    quad( 1, 0, 3, 2 );
    quad( 2, 3, 7, 6 );
    quad( 3, 0, 4, 7 );
    quad( 6, 5, 1, 2 );
    quad( 4, 5, 6, 7 );
    quad( 5, 4, 0, 1 );
}

//----------------------------------------------------------------------------

/* Define the three parts */
/* Note use of push/pop to return modelview matrix to its state before functions were entered and use
rotation, translation, and scaling to create instances of symbols (cube and cylinder) */

void base() {
    mat4 instance = ( Translate( 0.0, 0.5 * BASE_HEIGHT, 0.0 ) * Scale( BASE_WIDTH, BASE_HEIGHT, BASE_WIDTH ) );

    glUniformMatrix4fv( ModelView, 1, GL_TRUE, model_view * instance );
    glDrawArrays( GL_TRIANGLES, 0, NumRobot );
}

//----------------------------------------------------------------------------

void upper_arm() {
    mat4 instance = ( Translate( 0.0, 0.5 * UPPER_ARM_HEIGHT, 0.0 ) * Scale( UPPER_ARM_WIDTH, UPPER_ARM_HEIGHT, UPPER_ARM_WIDTH ) );

    glUniformMatrix4fv( ModelView, 1, GL_TRUE, model_view * instance );
    glDrawArrays( GL_TRIANGLES, 0, NumRobot );
}

//----------------------------------------------------------------------------

void lower_arm() {
    mat4 instance = ( Translate( 0.0, 0.5 * LOWER_ARM_HEIGHT, 0.0 ) * Scale( LOWER_ARM_WIDTH, LOWER_ARM_HEIGHT, LOWER_ARM_WIDTH ) );

    glUniformMatrix4fv( ModelView, 1, GL_TRUE, model_view * instance );
    glDrawArrays( GL_TRIANGLES, 0, NumRobot );
}

//----------------------------------------------------------------------------

void initSphere(){
  // Code in this function below are from Textbook 6th edition, page 60.
  float DegreesToRadians = M_PI / 180.0;

  // LineLoop
  int k = 0;
  for(float phi = -80.0; phi <= 80.0; phi += 20.0){
    float phir = phi * DegreesToRadians;
    float phir20 = (phi + 20.0) * DegreesToRadians;

    for(float theta = -180.0; theta <= 180.0; theta += 20.0){
      float thetar = theta * DegreesToRadians;
      sphereLoop_points[k] = point4(sin(thetar)*cos(phir), cos(thetar)*cos(phir), sin(phir), 1.0);
      k++;
      sphereLoop_points[k] = point4(sin(thetar)*cos(phir20), cos(thetar)*cos(phir20), sin(phir20), 1.0);
      k++;
    }
  }

  for(int i=0; i<NumSphereLoop; i++){
    sphereLoop_colors[i] = vertex_colors[0];
  }

  // Fan
  k = 0;

  sphereFan_points[k] = point4(0.0, 0.0, 1.0, 1.0);
  k++;

  float sin80 = sin(80.0*DegreesToRadians);
  float cos80 = cos(80.0*DegreesToRadians);

  for(float theta = -180; theta <= 180; theta += 20.0){
    float thetar = theta * DegreesToRadians;
    sphereFan_points[k] = point4(sin(thetar)*cos80, cos(thetar)*cos80, sin80, 1.0);
    k++;
  }

  sphereFan_points[k] = point4(0.0, 0.0, -1.0, 1.0);
  k++;

  for(float theta = -180; theta <= 180; theta += 20.0){
    float thetar = theta;
    sphereFan_points[k] = point4(sin(thetar)*cos80, cos(thetar)*cos80, sin80, 1.0);
    k++;
  }

  for(int i=0; i<NumSphereFan; i++){
    sphereFan_colors[i] = vertex_colors[0];
  }
}

//----------------------------------------------------------------------------

void getAngles(){
  // 3 angles for fetch
  fetchBaseAngle = (int)(360 - (atan2(old_z, old_x) * 180 / M_PI)) % 360;
  fetchBaseAngle = round(fetchBaseAngle);

  float disTo02 = sqrt(old_x*old_x + old_z*old_z + (old_y-2)*(old_y-2));
  float distance = sqrt(old_x*old_x + old_z*old_z + (old_y-7)*(old_y-7));
  if(distance == 5.25){
    // no need to move lower arm
    fetchLowerAngle = 0;
  } else{
    // if <0, LowerArm move anti-clockwise
    float angle1 = acos((5*5+disTo02*disTo02-5.25*5.25)/(2*5*disTo02)) * 180 / M_PI;
    float angle2 = asin(1.0*(old_y-2)/disTo02) * 180 / M_PI;
    fetchLowerAngle = 90-(angle1+angle2);
  }
  fetchLowerAngle = round(fetchLowerAngle);

  float angle3 = acos((5*5+5.25*5.25-disTo02*disTo02)/(2*5*5.25)) * 180 / M_PI;
  fetchUpperAngle = 180 - angle3;
  fetchUpperAngle = round(fetchUpperAngle);

  // 3 angles for put
  float newBaseAngle;
  float newLowerAngle;
  float newUpperAngle;
  newBaseAngle = (int)(360 - (atan2(new_z, new_x) * 180 / M_PI)) % 360;
  putBaseAngle = newBaseAngle - fetchBaseAngle; // if <0, use anti-clockwise
  putBaseAngle = round(putBaseAngle);

  disTo02 = sqrt(new_x*new_x + new_z*new_z + (new_y-2)*(new_y-2));
  distance = sqrt(new_x*new_x + new_z*new_z + (new_y-7)*(new_y-7));
  if(distance == 5.25){
    newLowerAngle = 0;
  } else{
    float angle1 = acos((5*5+disTo02*disTo02-5.25*5.25)/(2*5*disTo02)) * 180 / M_PI;
    float angle2 = asin(1.0*(new_y-2)/disTo02) * 180 / M_PI;
    newLowerAngle = 90-(angle1+angle2);
  }
  putLowerAngle = newLowerAngle - fetchLowerAngle; // if <0, use anti-clockwise
  putLowerAngle = round(putLowerAngle);

  angle3 = acos((5*5+5.25*5.25-disTo02*disTo02)/(2*5*5.25)) * 180 / M_PI;
  newUpperAngle = 180 - angle3;
  putUpperAngle = newUpperAngle - fetchUpperAngle; // if <0; use anti-clockwise
  putUpperAngle = round(putUpperAngle);
}



//----------------------------------------------------------------------------

void init( void ) {
  status = 1;

  // Load shaders and use the resulting shader program
  GLuint program = InitShader("vshader81.glsl", "fshader81.glsl");
  glUseProgram(program);

  GLuint vPosition = glGetAttribLocation(program, "vPosition");
  GLuint vColor = glGetAttribLocation(program, "vColor");

  ModelView = glGetUniformLocation(program, "ModelView");
  Projection = glGetUniformLocation(program, "Projection");

  glEnable(GL_DEPTH);
  glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);

  colorcube();

  // bind robot with vao and vbo
  glGenVertexArrays( 1, &vaos[0] );
  glBindVertexArray( vaos[0] );

  glGenBuffers( 1, &vbos[0] );
  glBindBuffer( GL_ARRAY_BUFFER, vbos[0] );

  glBufferData( GL_ARRAY_BUFFER, sizeof(robot_points) + sizeof(robot_colors), robot_points, GL_STATIC_DRAW );

  glBufferSubData( GL_ARRAY_BUFFER, 0, sizeof(robot_points), robot_points );
  glBufferSubData( GL_ARRAY_BUFFER, sizeof(robot_points), sizeof(robot_colors), robot_colors );

  glEnableVertexAttribArray( vPosition );
  glVertexAttribPointer( vPosition, 4, GL_FLOAT, GL_FALSE, 0, BUFFER_OFFSET(0) );

  glEnableVertexAttribArray( vColor );
  glVertexAttribPointer( vColor, 4, GL_FLOAT, GL_FALSE, 0, BUFFER_OFFSET(sizeof(robot_points)) );

  initSphere();
  getAngles();
  // bind sphere LineLoop with vao and vbo
  glGenVertexArrays( 1, &vaos[1] );
  glBindVertexArray( vaos[1] );

  glGenBuffers( 1, &vbos[1] );
  glBindBuffer( GL_ARRAY_BUFFER, vbos[1] );

  glBufferData( GL_ARRAY_BUFFER, sizeof(sphereLoop_points) + sizeof(sphereLoop_colors), sphereLoop_points, GL_STATIC_DRAW );

  glBufferSubData( GL_ARRAY_BUFFER, 0, sizeof(sphereLoop_points), sphereLoop_points );
  glBufferSubData( GL_ARRAY_BUFFER, sizeof(sphereLoop_points), sizeof(sphereLoop_colors), sphereLoop_colors );

  glEnableVertexAttribArray( vPosition );
  glVertexAttribPointer( vPosition, 4, GL_FLOAT, GL_FALSE, 0, BUFFER_OFFSET(0) );

  glEnableVertexAttribArray( vColor );
  glVertexAttribPointer( vColor, 4, GL_FLOAT, GL_FALSE, 0, BUFFER_OFFSET(sizeof(sphereLoop_points)) );

  // bind sphere Fan with vao and vbo
  glGenVertexArrays( 1, &vaos[2] );
  glBindVertexArray( vaos[2] );

  glGenBuffers( 1, &vbos[2] );
  glBindBuffer( GL_ARRAY_BUFFER, vbos[2] );

  glBufferData( GL_ARRAY_BUFFER, sizeof(sphereFan_points) + sizeof(sphereFan_colors), sphereFan_points, GL_STATIC_DRAW );

  glBufferSubData( GL_ARRAY_BUFFER, 0, sizeof(sphereFan_points), sphereFan_points );
  glBufferSubData( GL_ARRAY_BUFFER, sizeof(sphereFan_points), sizeof(sphereFan_colors), sphereFan_colors );

  glEnableVertexAttribArray( vPosition );
  glVertexAttribPointer( vPosition, 4, GL_FLOAT, GL_FALSE, 0, BUFFER_OFFSET(0) );

  glEnableVertexAttribArray( vColor );
  glVertexAttribPointer( vColor, 4, GL_FLOAT, GL_FALSE, 0, BUFFER_OFFSET(sizeof(sphereFan_points)) );

  glClearColor( 1.0, 1.0, 1.0, 1.0 );
}

//----------------------------------------------------------------------------

void display( void ) {
  glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );

  glBindVertexArray(vaos[0]);
  model_view = mat4(1.0);

  if(topview_flag){
    // move to top view
    model_view *= Translate(0.0, 5.0, 0.0);
    model_view *= RotateX(90);
  }

  model_view *= RotateY(Theta[Base]);
  base();

  model_view *= (Translate(0.0, BASE_HEIGHT, 0.0) * RotateZ(Theta[LowerArm]));
  lower_arm();

  model_view *= (Translate(0.0, LOWER_ARM_HEIGHT, 0.0) * RotateZ(Theta[UpperArm]));
  upper_arm();


  sphere_view = mat4(1.0);
  if(topview_flag){
    // move to top view
    sphere_view *= Translate(0.0, 5.0, 0.0);
    sphere_view *= RotateX(90);
  }

  if(status == 1)
    sphere_view *= Translate(old_x, old_y, old_z);

  if(status == 2){
    sphere_view = model_view;
    sphere_view *= Translate(0.0, 5.25, 0.0);
  }

  if(status == 3)
    sphere_view *= Translate(new_x, new_y, new_z);

  sphere_view *= Scale(UPPER_ARM_WIDTH/2, UPPER_ARM_WIDTH/2, UPPER_ARM_WIDTH/2);

  // draw sphere LineLoop with GL_LINE_LOOP
  glBindVertexArray(vaos[1]);
  glUniformMatrix4fv(ModelView, 1, GL_TRUE, sphere_view);
  glDrawArrays(GL_TRIANGLE_FAN, 0, NumSphereLoop);

  // draw sphere Fan with GL_TRIANGLE_FAN
  glBindVertexArray(vaos[2]);
  glUniformMatrix4fv(ModelView, 1, GL_TRUE, sphere_view);
  glDrawArrays(GL_TRIANGLE_FAN, 0, NumSphereFan);

  glutSwapBuffers();
}

//----------------------------------------------------------------------------

void armRotate(int direction) {
  // 0 means left(anti-clockwise), 1 means right(clockwise)
  if(direction == 0){
    // Incrase the joint angle
	  Theta[Axis] += 1.0;
    if(Theta[Axis] >= 360.0){
      Theta[Axis] -= 360.0;
    }
  }
  if(direction == 1){
    // Decrase the joint angle
	  Theta[Axis] -= 1.0;
    if(Theta[Axis] <= 0.0){
      Theta[Axis] += 360.0;
    }
  }
  glutPostRedisplay();
}

void fetchput(int) {
  if(status == 1){
    // fetch starts here
    if(basefetch_count < fetchBaseAngle){
      Axis = 0;
      armRotate(0);
      basefetch_count++;
    } else if(lowerfetch_count < abs(fetchLowerAngle)){
      Axis = 1;
      if(fetchLowerAngle<0){
        armRotate(0);
      } else{
        armRotate(1);
      }
      lowerfetch_count++;
    } else if(upperfetch_count < fetchUpperAngle){
      // here upperArm must rotate anti-clockwise
      Axis = 2;
      armRotate(1);
      upperfetch_count++;
    } else {
      // fetch stops here
      status = 2;
      printf("start putting the ball!\n");
    }
  } else if(status == 2){
    // put starts here
    if(baseput_count < abs(putBaseAngle)){
      Axis = 0;
      if(putBaseAngle<0){
        armRotate(1);
      } else{
        armRotate(0);
      }
      baseput_count++;
    } else if(lowerput_count < abs(putLowerAngle)){
      Axis = 1;
      if(putLowerAngle < 0){
        armRotate(0);
      } else{
        armRotate(1);
      }
      lowerput_count++;
    } else if(upperput_count < abs(putUpperAngle)){
      Axis = 2;
      if(putUpperAngle < 0){
        armRotate(0);
      } else{
        armRotate(1);
      }
      upperput_count++;
    } else{
      // put stops here
      status = 3;
      printf("start going back!\n");
    }
  } else if(status == 3){
    // go back starts here
    if(Theta[0] != 0){
      Axis = 0;
      armRotate(0);
    } else if(Theta[1] != 0){
      Axis = 1;
      armRotate(0);
    } else if(Theta[2] != 0){
      Axis = 2;
      armRotate(0);
    } else{
      // anime stops
      speed = -1;
      printf("stop moving!\n");
    }
  }

  glutTimerFunc(speed, fetchput, 0);
}

//----------------------------------------------------------------------------

void reshape( int width, int height ) {
    glViewport( 0, 0, width, height );

    GLfloat  left = -10.0, right = 10.0;
    GLfloat  bottom = -5.0, top = 15.0;
    GLfloat  zNear = -10.0, zFar = 10.0;

    GLfloat aspect = GLfloat(width)/height;

    if ( aspect > 1.0 ) {
      left *= aspect;
      right *= aspect;
    } else {
      bottom /= aspect;
      top /= aspect;
    }

    mat4 projection = Ortho(left, right, bottom, top, zNear, zFar);
    glUniformMatrix4fv(Projection, 1, GL_TRUE, projection);

    model_view = mat4(1.0);  // An Identity matrix
    sphere_view = mat4(1.0);
}

//----------------------------------------------------------------------------

void keyboard( unsigned char key, int x, int y ) {
  switch( key ) {
  case 033: // Escape Key
  case 'q':
  case 'Q':
    printf("Quit game!\n");
    exit( EXIT_SUCCESS );
	   break;
  case '1':
    printf("1 pressed, base selected!\n");
    Axis = Base;
    break;
  case '2':
    printf("2 pressed, lower arm selected!\n");
    Axis = LowerArm;
    break;
  case '3':
    printf("3 pressed, upper arm selected!\n");
    Axis = UpperArm;
    break;
  case 't':
    printf("t pressed, shift to top view\n");
    topview_flag = true;
    glutPostRedisplay();
    break;
  case 's':
    printf("s pressed, shift to side view\n");
    topview_flag = false;
    glutPostRedisplay();
    break;
  case 'v':
    printf("v pressed, change view\n");
    topview_flag = !topview_flag;
    glutPostRedisplay();
    break;
  }
}

void keyboardSpecial(int key, int x, int y) {
	switch (key) {
	case GLUT_KEY_LEFT:
		printf("left arrow pressed!\n");
		armRotate(0);
		break;
	case GLUT_KEY_RIGHT:
		printf("right arrow pressed!\n");
		armRotate(1);
		break;
	}
}

//----------------------------------------------------------------------------

int main(int argc, char **argv) {
    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);
    glutInitWindowSize(512, 512);
    glutInitContextVersion(3, 2);
    glutInitContextProfile(GLUT_CORE_PROFILE);
    glutCreateWindow("MyRobot!");

    // Iff you get a segmentation error at line 34, please uncomment the line below
    glewExperimental = GL_TRUE;
    glewInit();

    // parse command line argvs
    old_x = atof(argv[1]);
    old_y = atof(argv[2]);
    old_z = atof(argv[3]);
    new_x = atof(argv[4]);
    new_y = atof(argv[5]);
    new_z = atof(argv[6]);
    if(argc == 8){
      std::string view = argv[7];
      if(!view.compare("-tv"))
        topview_flag = true;
      if(!view.compare("-sv"))
        topview_flag = false;
    }
    init();

    glutDisplayFunc(display);
    glutReshapeFunc(reshape);
    glutTimerFunc(speed, fetchput, 0);
    glutKeyboardFunc(keyboard);
    glutSpecialFunc(keyboardSpecial);

    glutMainLoop();
    return 0;
}
