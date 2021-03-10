#include <iostream>
#include <math.h>
#include <random>
#include <vector>

#include "cloth.h"
#include "collision/plane.h"
#include "collision/sphere.h"

using namespace std;

Cloth::Cloth(double width, double height, int num_width_points,
             int num_height_points, float thickness) {
  this->width = width;
  this->height = height;
  this->num_width_points = num_width_points;
  this->num_height_points = num_height_points;
  this->thickness = thickness;

  buildGrid();
  buildClothMesh();
}

Cloth::~Cloth() {
  point_masses.clear();
  springs.clear();

  if (clothMesh) {
    delete clothMesh;
  }
}


bool isPinned(int x, int y, Cloth* cloth) {
    for (auto v : cloth->pinned) {
        if (v.at(0) == x && v.at(1) == y)
            return true;
    }
    return false;

}

double fRand(double fMin, double fMax) {
    double f = (double)rand() / RAND_MAX;
    return fMin + f * (fMax - fMin);
}

void Cloth::buildGrid() {
  // TODO (Part 1): Build a grid of masses and springs.
    double widthSpacing = width / num_width_points, heightSpacing = height / num_height_points;

    if (orientation == HORIZONTAL) {
        for (size_t i = 0; i < num_height_points; ++i) {
            for (size_t j = 0; j < num_width_points; ++j) {
                double x = i * widthSpacing, z = j * heightSpacing, y = 1.0;
                point_masses.emplace_back(PointMass(Vector3D(x, y, z), isPinned(i, j, this)));
            }
        }
    }

    else {
        srand(42);
        for (size_t i = 0; i < num_height_points; ++i) {
            for (size_t j = 0; j < num_width_points; ++j) {
                double x = i * widthSpacing,
                    z = fRand(-1.0 / 1000, 1.0 / 1000),
                    y = j * heightSpacing;
                point_masses.emplace_back(PointMass(Vector3D(x, y, z), isPinned(i, j, this)));
            }
        }
    }

    size_t vecSize = point_masses.size();
    PointMass* curr, * b;
    for (int i = 0; i < vecSize; ++i) {
        curr = &(point_masses.at(i));
        if (i - 1 >= 0 && (i / num_width_points == (i - 1) / num_width_points)) {
            b = &(point_masses.at(i - 1));
            springs.emplace_back(Spring(curr, b, e_spring_type::STRUCTURAL));
        }

        if (i - num_width_points >= 0) {
            b = &(point_masses.at(i - num_width_points));
            springs.emplace_back(Spring(curr, b, e_spring_type::STRUCTURAL));
        }
        
        if (i - 2 >= 0 && (i / num_width_points == (i - 2) / num_width_points)) {
            b = &(point_masses.at(i - 2));
            springs.emplace_back(Spring(curr, b, e_spring_type::BENDING));
        }

        if (i - 2 * num_width_points >= 0) {
            b = &(point_masses.at(i - 2 * num_width_points));
            springs.emplace_back(Spring(curr, b, e_spring_type::BENDING));
        }

        if (i % num_width_points != 0 && i - num_width_points - 1 >= 0) {
            b = &(point_masses.at(i - num_width_points - 1));
            springs.emplace_back(Spring(curr, b, e_spring_type::SHEARING));
        }

        if ((i + 1) % num_width_points != 0 && i - num_width_points + 1 >= 0) {
            b = &(point_masses.at(i - num_width_points + 1));
            springs.emplace_back(Spring(curr, b, e_spring_type::SHEARING));
        }
    }

}


void Cloth::simulate(double frames_per_sec, double simulation_steps, ClothParameters *cp,
                     vector<Vector3D> external_accelerations,
                     vector<CollisionObject *> *collision_objects) {
    double mass = width * height * cp->density / num_width_points / num_height_points;
    double delta_t = 1.0f / frames_per_sec / simulation_steps;

    // TODO (Part 2): Compute total force acting on each point mass.
    Vector3D totalExternalForce;
    for (Vector3D v : external_accelerations) {
        totalExternalForce += v;
    }
    totalExternalForce *= mass;

    for (PointMass &p : point_masses) {
        p.forces.x = 0.0;
        p.forces.y = 0.0;
        p.forces.z = 0.0;        
        p.forces += totalExternalForce;        
    }

    for (Spring &s : springs) {
        if (cp->enable_bending_constraints && s.spring_type == e_spring_type::BENDING) {
            //Vector3D fSpring = cp->ks * ((s.pm_a->position - s.pm_b->position).norm() - s.rest_length) * (s.pm_a->position - s.pm_b->position).unit();
            Vector3D fSpring = cp->ks * 0.2 * ((s.pm_a->position - s.pm_b->position).norm() - s.rest_length) * (s.pm_a->position - s.pm_b->position).unit();
            s.pm_a->forces -= fSpring;
            s.pm_b->forces += fSpring; // might have this backwards. swap signs if so.
        }
        if (cp->enable_shearing_constraints&& s.spring_type == e_spring_type::SHEARING) {
            Vector3D fSpring = cp->ks * ((s.pm_a->position - s.pm_b->position).norm() - s.rest_length) * (s.pm_a->position - s.pm_b->position).unit();
            s.pm_a->forces -= fSpring;
            s.pm_b->forces += fSpring; // might have this backwards. swap signs if so.
        }
        if (cp->enable_structural_constraints && s.spring_type == e_spring_type::STRUCTURAL) {
            Vector3D fSpring = cp->ks * ((s.pm_a->position - s.pm_b->position).norm() - s.rest_length) * (s.pm_a->position - s.pm_b->position).unit();
            s.pm_a->forces -= fSpring;
            s.pm_b->forces += fSpring; // might have this backwards. swap signs if so.
        }
    }
  // TODO (Part 2): Use Verlet integration to compute new point mass positions
    for (PointMass &p : point_masses) {
        if (!p.pinned) {
            Vector3D newPosition = p.position + (1.0 - cp->damping / 100.0) * (p.position - p.last_position) + (p.forces / mass) * pow(delta_t, 2.0); 
            p.last_position = p.position;
            p.position = newPosition;
        }
    }

    

  // TODO (Part 4): Handle self-collisions.
    build_spatial_map();
    for (PointMass& p : point_masses) {
        self_collide(p, simulation_steps);
    }

  // TODO (Part 3): Handle collisions with other primitives.
    for (CollisionObject *c : *collision_objects) {
        for (PointMass& p : point_masses) {
            c->collide(p);
        }
    }

  // TODO (Part 2): Constrain the changes to be such that the spring does not change
  // in length more than 10% per timestep [Provot 1995].

    for (Spring& s : springs) {
        double actualLength = (s.pm_b->position - s.pm_a->position).norm();
        if (actualLength > s.rest_length * 1.1 + 0.001) { // 110% of rest length plus small rounding buffer
            double adjustedHalfLength = (actualLength - 1.1 * s.rest_length) / 2.0;
            if (!s.pm_a->pinned && !s.pm_b->pinned) {
                Vector3D a = (s.pm_b->position - s.pm_a->position) * adjustedHalfLength;
                s.pm_a->position += a;
                s.pm_b->position -= a;
                actualLength = (s.pm_b->position - s.pm_a->position).norm();
                Vector3D p(actualLength);

            }

            else if (s.pm_a->pinned && !s.pm_b->pinned) {
                Vector3D a = (s.pm_b->position - s.pm_a->position) * adjustedHalfLength * 2;
                s.pm_b->position -= a;
            }

            else if (!s.pm_a->pinned && s.pm_b->pinned) {
                Vector3D a = (s.pm_b->position - s.pm_a->position) * adjustedHalfLength * 2;
                s.pm_a->position += a;
            }
        }
    }

}

void Cloth::build_spatial_map() {
      for (const auto &entry : map) {
        delete(entry.second);
      }
      map.clear();

    // TODO (Part 4): Build a spatial map out of all of the point masses.
      for (PointMass &p : point_masses) {
          float key = hash_position(p.position);
          vector<PointMass*> *v = map[key];
          if (v == nullptr) {
              v = new vector<PointMass*>;
          }
          v->push_back(&p);
          map[key] = v;
      }
}

void Cloth::self_collide(PointMass &pm, double simulation_steps) {
    // TODO (Part 4): Handle self-collision for a given point mass.
    float key = hash_position(pm.position);
    vector<PointMass*>* v = map[key];
    Vector3D correctionVector;
    size_t count = 0;


    for (PointMass* candidate : *v) {
        if (!(pm.position == candidate->position && pm.last_position == candidate->last_position && pm.forces == candidate->forces)) {
            double distance = (pm.position - candidate->position).norm();
            if (distance < 2 * thickness) {
                count++;
                correctionVector += (candidate->position - pm.position).unit() * (2 * thickness - distance);
            }
        }
    }

    if (count != 0) {
        correctionVector /= count;
        correctionVector /= simulation_steps;
        pm.position -= correctionVector;
    }
    
}

float Cloth::hash_position(Vector3D pos) {
  // TODO (Part 4): Hash a 3D position into a unique float identifier that represents membership in some 3D box volume.
    float w = 3.0 * width / num_width_points, h = 3.0 * height / num_height_points;
    float t = max(w, h);
    
    float x = pos.x / w;
    float y = pos.y / h;
    float z = pos.z / t;

    return (int) (x + y + z);

  
}

///////////////////////////////////////////////////////
/// YOU DO NOT NEED TO REFER TO ANY CODE BELOW THIS ///
///////////////////////////////////////////////////////

void Cloth::reset() {
  PointMass *pm = &point_masses[0];
  for (int i = 0; i < point_masses.size(); i++) {
    pm->position = pm->start_position;
    pm->last_position = pm->start_position;
    pm++;
  }
}

void Cloth::buildClothMesh() {
  if (point_masses.size() == 0) return;

  ClothMesh *clothMesh = new ClothMesh();
  vector<Triangle *> triangles;

  // Create vector of triangles
  for (int y = 0; y < num_height_points - 1; y++) {
    for (int x = 0; x < num_width_points - 1; x++) {
      PointMass *pm = &point_masses[y * num_width_points + x];
      // Get neighboring point masses:
      /*                      *
       * pm_A -------- pm_B   *
       *             /        *
       *  |         /   |     *
       *  |        /    |     *
       *  |       /     |     *
       *  |      /      |     *
       *  |     /       |     *
       *  |    /        |     *
       *      /               *
       * pm_C -------- pm_D   *
       *                      *
       */
      
      float u_min = x;
      u_min /= num_width_points - 1;
      float u_max = x + 1;
      u_max /= num_width_points - 1;
      float v_min = y;
      v_min /= num_height_points - 1;
      float v_max = y + 1;
      v_max /= num_height_points - 1;
      
      PointMass *pm_A = pm                       ;
      PointMass *pm_B = pm                    + 1;
      PointMass *pm_C = pm + num_width_points    ;
      PointMass *pm_D = pm + num_width_points + 1;
      
      Vector3D uv_A = Vector3D(u_min, v_min, 0);
      Vector3D uv_B = Vector3D(u_max, v_min, 0);
      Vector3D uv_C = Vector3D(u_min, v_max, 0);
      Vector3D uv_D = Vector3D(u_max, v_max, 0);
      
      
      // Both triangles defined by vertices in counter-clockwise orientation
      triangles.push_back(new Triangle(pm_A, pm_C, pm_B, 
                                       uv_A, uv_C, uv_B));
      triangles.push_back(new Triangle(pm_B, pm_C, pm_D, 
                                       uv_B, uv_C, uv_D));
    }
  }

  // For each triangle in row-order, create 3 edges and 3 internal halfedges
  for (int i = 0; i < triangles.size(); i++) {
    Triangle *t = triangles[i];

    // Allocate new halfedges on heap
    Halfedge *h1 = new Halfedge();
    Halfedge *h2 = new Halfedge();
    Halfedge *h3 = new Halfedge();

    // Allocate new edges on heap
    Edge *e1 = new Edge();
    Edge *e2 = new Edge();
    Edge *e3 = new Edge();

    // Assign a halfedge pointer to the triangle
    t->halfedge = h1;

    // Assign halfedge pointers to point masses
    t->pm1->halfedge = h1;
    t->pm2->halfedge = h2;
    t->pm3->halfedge = h3;

    // Update all halfedge pointers
    h1->edge = e1;
    h1->next = h2;
    h1->pm = t->pm1;
    h1->triangle = t;

    h2->edge = e2;
    h2->next = h3;
    h2->pm = t->pm2;
    h2->triangle = t;

    h3->edge = e3;
    h3->next = h1;
    h3->pm = t->pm3;
    h3->triangle = t;
  }

  // Go back through the cloth mesh and link triangles together using halfedge
  // twin pointers

  // Convenient variables for math
  int num_height_tris = (num_height_points - 1) * 2;
  int num_width_tris = (num_width_points - 1) * 2;

  bool topLeft = true;
  for (int i = 0; i < triangles.size(); i++) {
    Triangle *t = triangles[i];

    if (topLeft) {
      // Get left triangle, if it exists
      if (i % num_width_tris != 0) { // Not a left-most triangle
        Triangle *temp = triangles[i - 1];
        t->pm1->halfedge->twin = temp->pm3->halfedge;
      } else {
        t->pm1->halfedge->twin = nullptr;
      }

      // Get triangle above, if it exists
      if (i >= num_width_tris) { // Not a top-most triangle
        Triangle *temp = triangles[i - num_width_tris + 1];
        t->pm3->halfedge->twin = temp->pm2->halfedge;
      } else {
        t->pm3->halfedge->twin = nullptr;
      }

      // Get triangle to bottom right; guaranteed to exist
      Triangle *temp = triangles[i + 1];
      t->pm2->halfedge->twin = temp->pm1->halfedge;
    } else {
      // Get right triangle, if it exists
      if (i % num_width_tris != num_width_tris - 1) { // Not a right-most triangle
        Triangle *temp = triangles[i + 1];
        t->pm3->halfedge->twin = temp->pm1->halfedge;
      } else {
        t->pm3->halfedge->twin = nullptr;
      }

      // Get triangle below, if it exists
      if (i + num_width_tris - 1 < 1.0f * num_width_tris * num_height_tris / 2.0f) { // Not a bottom-most triangle
        Triangle *temp = triangles[i + num_width_tris - 1];
        t->pm2->halfedge->twin = temp->pm3->halfedge;
      } else {
        t->pm2->halfedge->twin = nullptr;
      }

      // Get triangle to top left; guaranteed to exist
      Triangle *temp = triangles[i - 1];
      t->pm1->halfedge->twin = temp->pm2->halfedge;
    }

    topLeft = !topLeft;
  }

  clothMesh->triangles = triangles;
  this->clothMesh = clothMesh;
}
