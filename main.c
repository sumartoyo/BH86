#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <signal.h>

#define NODE_EMPTY 0
#define NODE_INTERNAL 1
#define NODE_EXTERNAL 2

float G = 0.0001;
float THETA = 0.5;
float EPS_2 = 0.001;
float ETA = 0.1;

int n_bodies;
float *bodies_mass;
float *bodies_x;
float *bodies_y;
float *bodies_vx;
float *bodies_vy;

int n_nodes;
int *nodes_type;
int *nodes_body;
float *nodes_mass;
float *nodes_xmass;
float *nodes_ymass;
float *nodes_xmin;
float *nodes_ymin;
float *nodes_width;
int *nodes_nw;
int *nodes_ne;
int *nodes_sw;
int *nodes_se;

void init_bodies();
void init_nodes();
void free_bodies();
void free_nodes();

void create_body(int idx_body, float x, float y, float mass);
void create_node(float xmin, float ymin, float width);

void compute_force(int idx_body, int idx_node, float *fx, float *fy);
void compute_force_recursive(int idx_body, int idx_node, float *fx, float *fy, int children[]);
void update_pos(int idx_body, float fx, float fy);

void insert(int idx_node, int idx_body);
void insert_recursive(int idx_node, int idx_body);
void update_mass(int idx_node, int idx_body);
void branch(int idx_node);
void eliminate_empty(int idx_node);
void eliminate_empty_recursive(int idx_node, int children[]);
void make_tree();

void step(int epoch);

void print_tree(int idx_node, int tab, char *name) {
    int i, idx_body;
    if (idx_node > -1) {
        for (i = 0; i < tab; i++) {
            printf("\t");
        }
        printf("%s %s",
            name,
            nodes_type[idx_node] == NODE_INTERNAL ? "nin" : "nex"
        );
        idx_body = nodes_body[idx_node];
        if (idx_body > -1) {
            printf(" body x=%f y=%f m=%f\n",
                bodies_x[idx_body],
                bodies_y[idx_body],
                bodies_mass[idx_body]
            );
        } else {
            printf(" cx=%f cy=%f cm=%f\n",
                nodes_xmass[idx_node],
                nodes_ymass[idx_node],
                nodes_mass[idx_node]
            );
        }
        print_tree(nodes_nw[idx_node], tab+1, "nw");
        print_tree(nodes_ne[idx_node], tab+1, "ne");
        print_tree(nodes_sw[idx_node], tab+1, "sw");
        print_tree(nodes_se[idx_node], tab+1, "se");
    }
}

static volatile int keepRunning = 1;
void intHandler(int dummy) {
    keepRunning = 0;
}

int main(int argc, char *argv[]) {
    signal(SIGINT, intHandler);

    n_bodies = 2;
    init_bodies();

    create_body(0, -0.1, 0, 0.1);
    create_body(1, 0.1, 0, 0.1);

    // make_tree();
    // print_tree(0, 0, "root");

    step(0);

    free_bodies();

    return 0;
}

// implementations

void init_bodies() {
    long long int size_float = n_bodies*sizeof(float);
    bodies_mass = malloc(size_float);
    bodies_x = malloc(size_float);
    bodies_y = malloc(size_float);
    bodies_vx = malloc(size_float);
    bodies_vy = malloc(size_float);
}

void init_nodes() {
    long long int size_float, size_int;
    float sisa = n_bodies*1.0;
    float total = 0;

    while (sisa >= 1) {
        total += sisa;
        sisa = sisa/2;
    }
    n_nodes = (int) ceil(total);

    size_float = n_nodes*sizeof(float);
    size_int = n_nodes*sizeof(int);
    nodes_type = malloc(size_int);
    nodes_body = malloc(size_int);
    nodes_mass = malloc(size_float);
    nodes_xmass = malloc(size_float);
    nodes_ymass = malloc(size_float);
    nodes_xmin = malloc(size_float);
    nodes_ymin = malloc(size_float);
    nodes_width = malloc(size_float);
    nodes_nw = malloc(size_int);
    nodes_ne = malloc(size_int);
    nodes_sw = malloc(size_int);
    nodes_se = malloc(size_int);

    n_nodes = 0;
}

void free_bodies() {
    free(bodies_mass);
    free(bodies_x);
    free(bodies_y);
    free(bodies_vx);
    free(bodies_vy);
}

void free_nodes() {
    free(nodes_type);
    free(nodes_body);
    free(nodes_mass);
    free(nodes_xmass);
    free(nodes_ymass);
    free(nodes_xmin);
    free(nodes_ymin);
    free(nodes_width);
    free(nodes_nw);
    free(nodes_ne);
    free(nodes_sw);
    free(nodes_se);
}

void create_body(int idx_body, float x, float y, float mass) {
    bodies_x[idx_body] = x;
    bodies_y[idx_body] = y;
    bodies_mass[idx_body] = mass;
}

void create_node(float xmin, float ymin, float width) {
    nodes_type[n_nodes] = NODE_EXTERNAL;
    nodes_body[n_nodes] = -1;
    nodes_xmin[n_nodes] = xmin;
    nodes_ymin[n_nodes] = ymin;
    nodes_width[n_nodes] = width;
    nodes_nw[n_nodes] = -1;
    nodes_ne[n_nodes] = -1;
    nodes_sw[n_nodes] = -1;
    nodes_se[n_nodes] = -1;
    n_nodes += 1;
}

void compute_force(int idx_body, int idx_node, float *fx, float *fy) {
    float r, strength;
    float dx = nodes_xmass[idx_node]-bodies_x[idx_body];
    float dy = nodes_ymass[idx_node]-bodies_y[idx_body];

    *fx = 0;
    *fy = 0;

    if (fabs(dx) > 0 || fabs(dy) > 0) {
        r = sqrt(pow(dx, 2) + pow(dy, 2));
        strength = G*bodies_mass[idx_body]*nodes_mass[idx_node] /
            sqrt(pow(pow(r, 2)+EPS_2, 3));

        if (nodes_body[idx_node] > -1) {
            *fx = strength*dx;
            *fy = strength*dy;
        } else {
            if (nodes_width[idx_node]/r < THETA) {
                *fx = strength*dx;
                *fy = strength*dy;
            } else {
                compute_force_recursive(idx_body, idx_node, fx, fy, nodes_nw);
                compute_force_recursive(idx_body, idx_node, fx, fy, nodes_ne);
                compute_force_recursive(idx_body, idx_node, fx, fy, nodes_sw);
                compute_force_recursive(idx_body, idx_node, fx, fy, nodes_se);
            }
        }
    }
}

void compute_force_recursive(int idx_body, int idx_node, float *fx, float *fy, int children[]) {
    float fx_child, fy_child;
    int idx_child = children[idx_node];
    if (idx_child > -1) {
        compute_force(idx_body, idx_child, &fx_child, &fy_child);
        *fx += fx_child;
        *fy += fy_child;
    }
}

void update_pos(int idx_body, float fx, float fy) {
    float ax = fx/bodies_mass[idx_body],
          ay = fy/bodies_mass[idx_body],
          dvx = ax*ETA,
          dvy = ay*ETA;
    bodies_x[idx_body] += (bodies_vx[idx_body]*ETA)+(0.5*dvx*ETA);
    bodies_y[idx_body] += (bodies_vy[idx_body]*ETA)+(0.5*dvy*ETA);
    bodies_vx[idx_body] += dvx;
    bodies_vy[idx_body] += dvy;
}

void insert(int idx_node, int idx_body) {
    update_mass(idx_node, idx_body);
    if (nodes_type[idx_node] == NODE_INTERNAL) {
        insert_recursive(idx_node, idx_body);
    } else if (nodes_body[idx_node] > -1) {
        nodes_type[idx_node] = NODE_INTERNAL;
        branch(idx_node);
        insert_recursive(idx_node, nodes_body[idx_node]);
        insert_recursive(idx_node, idx_body);
        nodes_body[idx_node] = -1;
    } else {
        nodes_body[idx_node] = idx_body;
    }
}

void insert_recursive(int idx_node, int idx_body) {
    float width_half = nodes_width[idx_node]/2;
    char insert_ver = bodies_y[idx_body] < nodes_ymin[idx_node]+width_half ? 's' : 'n';
    char insert_hor = bodies_x[idx_body] < nodes_xmin[idx_node]+width_half ? 'w' : 'e';
    if (insert_ver == 'n') {
        if (insert_hor == 'w') {
            insert(nodes_nw[idx_node], idx_body);
        } else {
            insert(nodes_ne[idx_node], idx_body);
        }
    } else {
        if (insert_hor == 'w') {
            insert(nodes_sw[idx_node], idx_body);
        } else {
            insert(nodes_se[idx_node], idx_body);
        }
    }
}

void update_mass(int idx_node, int idx_body) {
    float weight_x = nodes_xmass[idx_node]*nodes_mass[idx_node],
          weight_y = nodes_ymass[idx_node]*nodes_mass[idx_node];
    weight_x += bodies_x[idx_body]*bodies_mass[idx_body];
    weight_y += bodies_y[idx_body]*bodies_mass[idx_body];
    nodes_mass[idx_node] += bodies_mass[idx_body];
    nodes_xmass[idx_node] = weight_x/nodes_mass[idx_node];
    nodes_ymass[idx_node] = weight_y/nodes_mass[idx_node];
}

void branch(int idx_node) {
    float xmin = nodes_xmin[idx_node];
    float ymin = nodes_ymin[idx_node];
    float width_half = nodes_width[idx_node]/2;
    create_node(xmin, ymin+width_half, width_half);
    nodes_nw[idx_node] = n_nodes-1;
    create_node(xmin+width_half, ymin+width_half, width_half);
    nodes_ne[idx_node] = n_nodes-1;
    create_node(xmin, ymin, width_half);
    nodes_sw[idx_node] = n_nodes-1;
    create_node(xmin+width_half, ymin, width_half);
    nodes_se[idx_node] = n_nodes-1;
}

void eliminate_empty(int idx_node) {
    eliminate_empty_recursive(idx_node, nodes_nw);
    eliminate_empty_recursive(idx_node, nodes_ne);
    eliminate_empty_recursive(idx_node, nodes_sw);
    eliminate_empty_recursive(idx_node, nodes_se);
}

void eliminate_empty_recursive(int idx_node, int children[]) {
    int idx_child = children[idx_node];
    if (idx_child > -1) {
        if (nodes_type[idx_child] == NODE_EXTERNAL) {
            if (nodes_body[idx_child] == -1) {
                children[idx_node] = -1;
            }
        } else {
            eliminate_empty(idx_child);
        }
    }
}

void make_tree() {
    int i;
    float xmin = bodies_x[0],
          xmax = bodies_x[0],
          ymin = bodies_y[0],
          ymax = bodies_y[0];
    for (i = 1; i < n_bodies; i++) {
        if (bodies_x[i] < xmin) xmin = bodies_x[i];
        if (bodies_x[i] > xmax) xmax = bodies_x[i];
        if (bodies_y[i] < ymin) ymin = bodies_y[i];
        if (bodies_y[i] > ymax) ymax = bodies_y[i];
    }

    init_nodes();
    create_node(xmin, ymin, xmax-xmin > ymax-ymin ? xmax-xmin : ymax-ymin);
    for (i = 0; i < n_bodies; i++) {
        insert(0, i);
    }
    eliminate_empty(0);
}

void step(int epoch) {
    int i;
    float fx, fy;
    printf("epoch %d\n", epoch);
    make_tree();
    for (i = 0; i < n_bodies; i++) {
        compute_force(i, 0, &fx, &fy);
        update_pos(i, fx, fy);
        printf("%f %f %f\n", bodies_mass[i], bodies_x[i], bodies_vx[i]);
    }
    if (keepRunning) {
        usleep(50000);
        step(epoch+1);
    }
    free_nodes();
}
