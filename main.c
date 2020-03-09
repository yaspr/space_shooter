/*
  
 */

#include <math.h>
#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <sys/time.h>

#include "flame.h"

#include "pi.h"

#define max(a, b) ((a) > (b)) ? (a) : (b)

#define COLOR_WHITE 0xFF
#define COLOR_BLACK 0x00

#define MAX_STARS 300
#define MAX_BOMBS 300
#define MAX_LIVES 20

#define MIN_X 20
#define MAX_X 1800

#define MIN_Y 20
#define MAX_Y 1000

#define MAX_FLEET 200

#define MAX_BULLETS 30

#define MAX_DIST 300

#define MAX_LEVEL 100

#define MAX_STR 128

//
typedef unsigned char byte;

//
typedef struct particle_s {

  byte state;
  
  //Position
  float x;
  float y;
  
  //
  float vx;
  float vy;

  //
  int dist;
  int max_dist;
  
  //
  byte R, G, B;
  
} particle_t;

//
typedef struct vehicle_s {

  //
  byte state;
  
  //Position
  float x;
  float y;

  //
  float r;
  
  //
  float a;

  //
  float vx;
  float vy;

  //
  float max_force;
  float max_velocity;

  //
  byte refueling;
  
  //
  byte shooting;
  
  //
  byte shield;       //Active (1) or Down (0)
  float  shield_level; //1 .. 100

  //Shield color
  byte s_R, s_G, s_B;
  
  //
  int nb_bullets;
  int max_bullets;

  particle_t *bullets;

  //
  int nb_kills;

  //
  double s_step;
  
  //
  byte R, G, B;

} vehicle_t;


//
static inline unsigned randxy(unsigned x, unsigned y)
{ return (rand() % (y - x + 1)) + x; }

//Distance between two points (Pythagoras)
static inline float dist(float x1, float y1, float x2, float y2)
{
  return sqrt((x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2));
}

//
static inline void norm(float *x, float *y)
{
  float d = sqrt((*x) * (*x) + (*y) * (*y));

  (*x) /= d;
  (*y) /= d;
}

//
static inline void scale(float *x, float *y, float s)
{
  *x *= s;
  *y *= s;
}

//
static inline void limit(float *x, float *y, float s)
{
  norm(x, y);
  scale(x, y, s);
}

//
void draw_vehicle(flame_obj_t *fo, vehicle_t *v, byte color)
{
  //
  if (color)
    flame_set_color(fo, v->R, v->G, v->B);
  else
    flame_set_color(fo, 0x00, 0x00, 0x00);

  //
  flame_draw_point(fo, v->x, v->y);
  
  //Triangle direction
  flame_draw_line(fo, v->x, v->y, v->x + v->r * cos(v->a), v->y + v->r * sin(v->a));

  if (v->shield)
    {
      if (color)
	flame_set_color(fo, v->s_R, v->s_G, v->s_B);
      else
	flame_set_color(fo, 0x00, 0x00, 0x00);
      
      for (double a = 0.0; a < 2 * PI; a += v->s_step)
	flame_draw_point(fo, v->x + 2 * v->r * cos(a), v->y + 2 * v->r * sin(a));
    }
  
  //
  if (color)
    flame_set_color(fo, v->R, v->G, v->B);
  else
    flame_set_color(fo, 0x00, 0x00, 0x00);
  
  //Delta vehicle (Triangle)
  flame_draw_line(fo,
		  v->x + v->r * cos(v->a),
		  v->y + v->r * sin(v->a),
		  v->x + v->r * cos(PI_2 + PI_3 + v->a),
		  v->y + v->r * sin(PI_2 + PI_3 + v->a));

  flame_draw_line(fo,
		  v->x + v->r * cos(v->a),
		  v->y + v->r * sin(v->a),
		  v->x + v->r * cos(PI + PI_6 + v->a),
		  v->y + v->r * sin(PI + PI_6 + v->a));

  flame_draw_line(fo,
		  v->x + v->r * cos(PI_2 + PI_3 + v->a),
		  v->y + v->r * sin(PI_2 + PI_3 + v->a),
		  v->x + v->r * cos(PI + PI_6 + v->a),
		  v->y + v->r * sin(PI + PI_6 + v->a));
}

//
void draw_particle(flame_obj_t *fo, particle_t *p, byte color)
{
  //
  if (color)
    flame_set_color(fo, p->R, p->G, p->B);
  else
    flame_set_color(fo, 0x00, 0x00, 0x00);
  
  flame_draw_point(fo, p->x, p->y);  
}

//
void update_vehicle(vehicle_t *v, float dt)
{
  v->x += v->vx * dt; 
  v->y += v->vy * dt;

  /* if (v->x <= MIN_X || v->x >= MAX_X) */
  /*   v->vx *= -1; */
  
  /* if (v->y <= MIN_Y || v->y >= MAX_Y) */
  /*   v->vy *= -1; */
}  

//
void update_particle(particle_t *p, float dt)
{
  p->x += p->vx * dt;
  p->y += p->vy * dt;
}

//
void seek_vehicle(vehicle_t *v, float tx, float ty, float dt)
{
  float dist_x = tx - v->x;
  float dist_y = ty - v->y;

  limit(&dist_x, &dist_y, v->max_velocity);
  
  float steer_x = dist_x - (v->vx * dt);
  float steer_y = dist_y - (v->vy * dt);

  limit(&steer_x, &steer_y, v->max_force);
  
  v->vx += steer_x;
  v->vy += steer_y;
  
  limit(&v->vx, &v->vy, v->max_velocity);
  
  if (v->y < ty)
    v->a = PI_2 - atan(dist_x / dist_y);
  else
    v->a = -PI_2 - atan(dist_x / dist_y);  
}

//
int main(int argc, char **argv)
{
  char c;

  //Event for keyboard & mouse interruptions
  XEvent event;

  //Set display size
  int x_min = 0, x_max = 1920;
  int y_min = 0, y_max = 1080;
  
  //Initialize display
  flame_obj_t *fo = flame_open("Delta", x_max, y_max);

  //
  int scale = 1;
  int bx = x_max >> 1, by = y_max >> 1; //Start in the middle of the screen

  //
  double dt = 0.0;
  double pdt = 0.0;
  double elapsed = 0.0;
  double time_count = 0.0;
  double time_step_s = 0.0;
  double _before_, _after_;
  struct timeval before, after;
  unsigned long long frame_count = 0;

  //
  char header[MAX_STR];
  
  //
  srand(getpid());

  //
  vehicle_t v;
  vehicle_t m;
  vehicle_t fleet[MAX_FLEET];
  
 lbl_start:
  
  //UFO
  v.x = 600;
  v.y = 600;
  v.r = 30;
  v.a = -PI_2;

  v.vx = 1.0;
  v.vy = 1.0;
  
  v.shield = 1;
  v.shield_level = MAX_LEVEL;
  
  v.nb_bullets = 0;
  v.max_bullets = MAX_BULLETS * 10;

  v.bullets = malloc(sizeof(particle_t) * MAX_BULLETS * 10);

  v.nb_kills = 0;

  v.s_step = 0.01;
  
  v.s_R = 0;
  v.s_G = 0xFF;
  v.s_B = 0;
  
  v.R = v.G = v.B = 0xFF;

  //Mothership
  m.x = randxy(MIN_X, MAX_X);
  m.y = randxy(MIN_Y, MAX_Y);

  m.r = 60;
  
  m.a = 0;

  m.vx = 0.5;
  m.vy = 0.5;
  
  m.shield = 1;
  
  m.R = 0xFF;
  m.G = 0;
  m.B = 0;

  m.s_R = 0xFF;
  m.s_G = 0;
  m.s_B = 0;

  m.s_step = 0.01;
  
  //Drones fleet
  fleet[0].state = 1;

  fleet[0].x = m.x;
  fleet[0].y = m.y;
  fleet[0].r = 10;
  
  fleet[0].a = 0.0;
      
  fleet[0].vx = 1;
  fleet[0].vy = 1;

  fleet[0].max_force = 0.8;
  fleet[0].max_velocity = 10;

  fleet[0].shield = 0;

  fleet[0].refueling = 0;
  
  fleet[0].nb_bullets = 0;
  fleet[0].max_bullets = MAX_BULLETS;

  fleet[0].bullets = malloc(sizeof(particle_t));
  
  fleet[0].bullets[0].dist = 0;
  fleet[0].bullets[0].max_dist = MAX_DIST * 3;

  fleet[0].R = 0;
  fleet[0].G = 0;
  fleet[0].B = 0xFF;

  double a = 0.0;
  
  for (unsigned i = 1; i < MAX_FLEET; i++, a += 1.0)
    {
      fleet[i].state = 1;
      
      fleet[i].x = m.x + randxy(2, 5) * m.r * cos(a);
      fleet[i].y = m.y + randxy(2, 5) * m.r * sin(a);
      fleet[i].r = 10;
      fleet[i].a = 0.0;
      
      fleet[i].vx = 1;
      fleet[i].vy = 1;
      
      fleet[i].max_force = 0.8;
      fleet[i].max_velocity = 10;

      fleet[i].shield = 0;

      fleet[i].refueling = 0;
      
      fleet[i].nb_bullets = 0;
      fleet[i].max_bullets = MAX_BULLETS;
      
      fleet[i].bullets = malloc(sizeof(particle_t));
      
      fleet[i].bullets[0].dist = 0;
      fleet[i].bullets[0].max_dist = MAX_DIST * 3;
      
      fleet[i].R = 0;
      fleet[i].G = 0;
      fleet[i].B = 0xFF;
    }
  
  //
  flame_clear_display(fo);

  //
  gettimeofday(&before, NULL);
  
  //
  while (1)
    {
      gettimeofday(&after, NULL);

      //Time in ms
      _before_ = (before.tv_sec * 1000) + (before.tv_usec) / 1000.0;
      _after_  = (after.tv_sec  * 1000) + (after.tv_usec)  / 1000.0;

      elapsed = (_after_ - _before_);
      
      time_count += elapsed;

      dt = (elapsed / 200.0);
      
      before   = after;
      _before_ = _after_;

      frame_count++;
      
      //Handle input events
      if (XPending(fo->display) > 0)
	{
	  XNextEvent(fo->display, &event);

	  //Keyboard input
	  if (event.type == KeyPress)
	    {
	      c = XLookupKeysym(&event.xkey, 0);
	      
	      if (c == 'q')
		break;
	      else
		if (c == 'r')
		  {
		    free(v.bullets);

		    for (unsigned i = 0; i < MAX_FLEET; i++)
		      free(fleet[i].bullets);
		    
		    goto lbl_start;
		  }
		else
		  {
		    switch (c)
		      {			
			//Left
		      case 81:
			draw_vehicle(fo, &v, 0);
			
			v.a -= 0.08;
			
			draw_vehicle(fo, &v, 1);
			break;

			//Right
		      case 83:
			draw_vehicle(fo, &v, 0);

			v.a += 0.08;
			
			draw_vehicle(fo, &v, 1);
			break;

			//Up
		      case 82:
			draw_vehicle(fo, &v, 0);

			if (v.x > MIN_X && v.x < MAX_X && v.y > MIN_Y && v.y < MAX_Y)
			  {
			    v.x += (v.r * cos(v.a)) / 2;
			    v.y += (v.r * sin(v.a)) / 2;
			  }
			else
			  {
			    v.x -= (v.r * cos(v.a)) / 2;
			    v.y -= (v.r * sin(v.a)) / 2;
			  }
			
			draw_vehicle(fo, &v, 1);
			break;

			//Down
		      case 84:
			draw_vehicle(fo, &v, 0);

			if (v.x > MIN_X && v.x < MAX_X && v.y > MIN_Y && v.y < MAX_Y)
			  {
			    v.x -= (v.r * cos(v.a)) / 2;
			    v.y -= (v.r * sin(v.a)) / 2;
			  }
			else
			  {
			    v.x += (v.r * cos(v.a)) / 2;
			    v.y += (v.r * sin(v.a)) / 2;
			  }
			
			draw_vehicle(fo, &v, 1);
			break;

			//Renew shield
		      case 's':
			v.shield = 1;
			v.shield_level = 100.0;

			v.s_R = 0;
			v.s_G = 0xFF;
			v.s_B = 0;
			
			break;
			
			//Shooting
		      case ' ':
			
			if (v.nb_bullets < v.max_bullets)
			  {
			    v.shooting = 1;

			    v.bullets[v.nb_bullets].state = 1;
			    
			    v.bullets[v.nb_bullets].x = v.x;
			    v.bullets[v.nb_bullets].y = v.y;
			    
			    v.bullets[v.nb_bullets].vx = v.r * cos(v.a) * 5;
			    v.bullets[v.nb_bullets].vy = v.r * sin(v.a) * 5;

			    v.bullets[v.nb_bullets].dist = 0;
			    v.bullets[v.nb_bullets].max_dist = MAX_DIST / dt;
			    
			    v.bullets[v.nb_bullets].R = 0xFF;
			    v.bullets[v.nb_bullets].G = 0xFF;
			    v.bullets[v.nb_bullets].B = 0xFF;
			    
			    v.nb_bullets++;
			  }
			break;
		      }
		  }
	    }
	  else //Mouse input
	    if (event.type == ButtonPress)
	      {
	      }
	}

      //Main game loop
      
      //Draw spacecraft
      draw_vehicle(fo, &v, 0);
      update_vehicle(&v, dt);
      draw_vehicle(fo, &v, 1);

      draw_vehicle(fo, &m, 0);
      update_vehicle(&m, dt);
      draw_vehicle(fo, &m, 1);

      //
      if (v.nb_kills == MAX_FLEET)
	{
	  printf("All drones were shut down\n");

	  goto lbl_start;
	}
      
      //
      if (v.shooting)
	{
	  for (unsigned i = 0; i < v.nb_bullets; i++)
	    {
	      if (v.bullets[i].state)
	      	{
		  if (v.bullets[i].dist < v.bullets[i].max_dist)
		    {
		      draw_particle(fo, &v.bullets[i], 0);
		      update_particle(&v.bullets[i], dt);
		      draw_particle(fo, &v.bullets[i], 1);
		      
		      v.bullets[i].dist++;

		      //Check collision with an attack drone
		      //If drone not killed and if not refueling
		      for (unsigned j = 0; j < MAX_FLEET; j++)
			if (fleet[j].state && !fleet[j].refueling)
			  if (dist(v.bullets[i].x, v.bullets[i].y, fleet[j].x, fleet[j].y) < fleet[j].r)
			    {
			      //Drone is shot
			      fleet[j].state = 0;
			      
			      draw_vehicle(fo, &fleet[j], 0);
			      
			      //Bullet is done
			      v.bullets[i].state = 0;
			      
			      draw_particle(fo, &v.bullets[i], 0);

			      v.nb_kills++;
			    }
		    }
		  else
		    {
		      v.bullets[i].state = 0;
		      
		      draw_particle(fo, &v.bullets[i], 0);
		    }
		}
	      else
		{
		  draw_particle(fo, &v.bullets[i], 0);	
		}
	    }
	}
      
      //Handle drone fleet
      for (unsigned i = 0; i < MAX_FLEET; i++)
	{
	  if (fleet[i].state)
	    {
	      draw_vehicle(fo, &fleet[i], 0);

	      if (fleet[i].nb_bullets >= fleet[i].max_bullets)
		{
		  fleet[i].refueling = 1;

		  fleet[i].R = 0;
		  fleet[i].G = 0xFF;
		  fleet[i].B = 0;
		  
		  if (dist(fleet[i].x, fleet[i].y, m.x, m.y) > (2 * m.r + fleet[i].r))
		    seek_vehicle(&fleet[i], m.x, m.y, dt);
		  else
		    {
		      //Drone fueling time
		      if (time_count >= 2000)
			{
			  fleet[i].refueling = 0;
			  fleet[i].nb_bullets = 0;
			  
			  seek_vehicle(&fleet[i], v.x, v.y, dt);

			  fleet[i].R = 0;
			  fleet[i].G = 0;
			  fleet[i].B = 0xFF;

			  time_count = 0;
			}
		      else
			{
			  fleet[i].vx += 2 * m.r * cos(m.a) * dt;
			  fleet[i].vy += 2 * m.r * sin(m.a) * dt;
			}
		    }
		}
	      
	      //If distant from target
	      if (dist(fleet[i].x, fleet[i].y, v.x, v.y) > (10 * v.r + fleet[i].r))
		{
		  if (fleet[i].nb_bullets < fleet[i].max_bullets)
		    seek_vehicle(&fleet[i], v.x, v.y, dt);

		  if (!fleet[i].refueling)
		    {
		      fleet[i].R = 0;
		      fleet[i].G = 0;
		      fleet[i].B = 0xFF;
		    }
		}
	      else
		{
		  //
		  if (!fleet[i].bullets[0].state)
		    {
		      //
		      if (fleet[i].nb_bullets < fleet[i].max_bullets)
			{		  
			  fleet[i].R = 0xFF;
			  fleet[i].G = 0;
			  fleet[i].B = 0;
		      
			  fleet[i].bullets[0].x = fleet[i].x;
			  fleet[i].bullets[0].y = fleet[i].y;
			  
			  fleet[i].bullets[0].vx = (v.x - fleet[i].x); 
			  fleet[i].bullets[0].vy = (v.y - fleet[i].y);
		      
			  fleet[i].bullets[0].dist = 0;

			  fleet[i].bullets[0].state = 1;
			  
			  fleet[i].bullets[0].R = 0xFF;
			  fleet[i].bullets[0].G = 0xFF;
			  fleet[i].bullets[0].B = 0xFF;
			  
			  fleet[i].nb_bullets++;
			}
		    }
		}
	      
	      //
	      if (fleet[i].bullets[0].state)
		{
		  //Make sure the bullet runs for a certain distance
		  if (fleet[i].bullets[0].dist < fleet[i].bullets[0].max_dist)
		    {
		      draw_particle(fo, &fleet[i].bullets[0], 0);
		      update_particle(&fleet[i].bullets[0], dt);
		      draw_particle(fo, &fleet[i].bullets[0], 1);
		  
		      //If bullet collides with target
		      if (dist(fleet[i].bullets[0].x, fleet[i].bullets[0].y, v.x, v.y) <= v.r)
			{
			  draw_particle(fo, &fleet[i].bullets[0], 0);

			  fleet[i].bullets[0].state = 0;
			  
			  //
			  if (v.shield_level > 0.0)
			    v.shield_level -= (2.0 * dt);

			  //Less than 50%
			  if (v.shield_level > 0.0 && v.shield_level < (MAX_LEVEL / 2.0))
			    {
			      //Turn the shield red
			      v.s_R = 0xFF;
			      v.s_G = 0;
			      v.s_B = 0;
			    }
		      
			  //
			  if (v.shield_level <= 0)
			    {
			      //Drop shield
			      v.shield_level = 0.0;
				  
			      draw_vehicle(fo, &v, 0);

			      v.shield = 0;
			    }		      
			}

		      fleet[i].bullets[0].dist++;
		    }
		  else
		    {
		      //
		      fleet[i].bullets[0].state = 0;
		      
		      draw_particle(fo, &fleet[i].bullets[0], 0);
		    }
		}
	      
	      update_vehicle(&fleet[i], dt);
	      draw_vehicle(fo, &fleet[i], 1);
	    }
	}
      
      //
      sprintf(header, "Shield: %10.3f    Bullets: %5d / %5d    Nb kills: %d", v.shield_level, v.nb_bullets, v.max_bullets, v.nb_kills);
      
      flame_set_color(fo, 0xFF, 0xFF, 0xFF);
      flame_draw_text(fo, 100, 100, header);
      
    }
  
  //
  for (unsigned i = 0; i < MAX_FLEET; i++)
    free(fleet[i].bullets);
  
  free(v.bullets);
  
  //
  flame_close(fo);

  return 0;
}
