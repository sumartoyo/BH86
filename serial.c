/* nama: Dimas Gilang
 * nim:  23515027
 * penggunaan:
 *     kompilasi:     gcc serial.c -o serial -lm
 *     melihat usage: ./serial
 *     contoh:        ./serial -n_bodies 1024
 * keterangan:
 *     kode ini tidak well commented
 *     penjelasan kode dapat dibaca di laporan
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <signal.h>

#define true 1
#define false 0

#define NODE_EMPTY 0
#define NODE_INTERNAL 1
#define NODE_EXTERNAL 2

// parameters
float G = 0.001;
float THETA = 0.1;
float EPS2 = 0.001;
float DTIME = 0.1;

// body data structure
int n_bodies;
float *bodies_mass;
float *bodies_x;
float *bodies_y;
float *bodies_vx;
float *bodies_vy;

// node data structure
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

// program variables
static volatile int steps = 10;
int is_print_bodies = false;
int is_print_tree = false;

// body functions
void create_body_random();
void create_body(int idx_body, float x, float y, float mass);
void compute_force(int idx_body, int idx_node, float *fx, float *fy);
void compute_force_recursive(int idx_body, int idx_node, float *fx, float *fy, int children[]);
void update_pos(int idx_body, float fx, float fy);
void print_bodies(char *title);

// node functions
void create_node(float xmin, float ymin, float width);
void insert(int idx_node, int idx_body);
void insert_recursive(int idx_node, int idx_body);
void update_mass(int idx_node, int idx_body);
void branch(int idx_node);
void eliminate_empty(int idx_node);
void eliminate_empty_recursive(int idx_node, int children[]);
void make_tree();
void print_tree(int idx_node, int tab, char *name);

// main loop
void step();
void handle_sigint(int _);

// program
int read_arg(int argc, char *argv[]);
int read_data();

// malloc/free
void malloc_bodies();
void malloc_nodes();
void free_bodies();
void free_nodes();

// main function
int main(int argc, char *argv[]) {
    if (read_arg(argc, argv) == false) {
        return 1;
    }
    signal(SIGINT, handle_sigint);
    if (is_print_bodies) print_bodies("[INPUT]");
    step();
    free_bodies();
    return 0;
}

// implementations

// body functions

void create_body_random() {
    int i;
    float x, y, mass;

    srand(time(NULL));
    for (i = 0; i < n_bodies; i++) {
        x = rand()%n_bodies;
        y = rand()%n_bodies;
        x += (rand()%10)*0.1;
        y += (rand()%10)*0.1;
        x *= rand()%2 == 0 ? 1 : -1;
        y *= rand()%2 == 0 ? 1 : -1;
        mass = (rand()%9)+11;
        create_body(i, x, y, mass);
    }
}

void create_body(int i, float x, float y, float mass) {
    bodies_x[i] = x;
    bodies_y[i] = y;
    bodies_mass[i] = mass;
}

void compute_force(int i, int j, float *fx, float *fy) {
    float r, strength;
    float dx = nodes_xmass[j]-bodies_x[i];
    float dy = nodes_ymass[j]-bodies_y[i];

    *fx = 0;
    *fy = 0;

    if (fabs(dx) > 0 || fabs(dy) > 0) {
        r = sqrt(pow(dx, 2) + pow(dy, 2));
        strength = G*bodies_mass[i]*nodes_mass[j] /
            sqrt(pow(pow(r, 2)+EPS2, 3));

        if (nodes_body[j] > -1) {
            *fx = strength*dx;
            *fy = strength*dy;
        } else {
            if (nodes_width[j]/r < THETA) {
                *fx = strength*dx;
                *fy = strength*dy;
            } else {
                compute_force_recursive(i, j, fx, fy, nodes_nw);
                compute_force_recursive(i, j, fx, fy, nodes_ne);
                compute_force_recursive(i, j, fx, fy, nodes_sw);
                compute_force_recursive(i, j, fx, fy, nodes_se);
            }
        }
    }
}

void compute_force_recursive(int i, int j, float *fx, float *fy, int children[]) {
    float fx_child, fy_child;
    int idx_child = children[j];
    if (idx_child > -1) {
        compute_force(i, idx_child, &fx_child, &fy_child);
        *fx += fx_child;
        *fy += fy_child;
    }
}

void update_pos(int i, float fx, float fy) {
    float ax = fx/bodies_mass[i],
          ay = fy/bodies_mass[i],
          dvx = ax*DTIME,
          dvy = ay*DTIME;
    bodies_vx[i] += dvx;
    bodies_vy[i] += dvy;
    bodies_x[i] += (bodies_vx[i]*DTIME)+(0.5*dvx*DTIME);
    bodies_y[i] += (bodies_vy[i]*DTIME)+(0.5*dvy*DTIME);
}

void print_bodies(char *title) {
    int i;
    printf("%s\n", title);
    printf("x  y  mass  vx  vy\n");
    for (i = 0; i < n_bodies; i++) {
        printf("%f %f %f %f %f\n",
            bodies_x[i],
            bodies_y[i],
            bodies_mass[i],
            bodies_vx[i],
            bodies_vy[i]
        );
    }
}

// node functions

void create_node(float xmin, float ymin, float width) {
    int j = n_nodes;
    nodes_type[j] = NODE_EXTERNAL;
    nodes_body[j] = -1;
    nodes_xmin[j] = xmin;
    nodes_ymin[j] = ymin;
    nodes_width[j] = width;
    nodes_nw[j] = -1;
    nodes_ne[j] = -1;
    nodes_sw[j] = -1;
    nodes_se[j] = -1;
    n_nodes += 1;
}

void insert(int j, int i) {
    update_mass(j, i);
    if (nodes_type[j] == NODE_INTERNAL) {
        insert_recursive(j, i);
    } else if (nodes_body[j] > -1) {
        nodes_type[j] = NODE_INTERNAL;
        branch(j);
        insert_recursive(j, nodes_body[j]);
        insert_recursive(j, i);
        nodes_body[j] = -1;
    } else {
        nodes_body[j] = i;
    }
}

void insert_recursive(int j, int i) {
    float width_half = nodes_width[j]/2;
    char insert_ver = bodies_y[i] < nodes_ymin[j]+width_half ? 's' : 'n';
    char insert_hor = bodies_x[i] < nodes_xmin[j]+width_half ? 'w' : 'e';
    if (insert_ver == 'n') {
        if (insert_hor == 'w') {
            insert(nodes_nw[j], i);
        } else {
            insert(nodes_ne[j], i);
        }
    } else {
        if (insert_hor == 'w') {
            insert(nodes_sw[j], i);
        } else {
            insert(nodes_se[j], i);
        }
    }
}

void update_mass(int j, int i) {
    float weight_x = nodes_xmass[j]*nodes_mass[j],
          weight_y = nodes_ymass[j]*nodes_mass[j];
    weight_x += bodies_x[i]*bodies_mass[i];
    weight_y += bodies_y[i]*bodies_mass[i];
    nodes_mass[j] += bodies_mass[i];
    nodes_xmass[j] = weight_x/nodes_mass[j];
    nodes_ymass[j] = weight_y/nodes_mass[j];
}

void branch(int j) {
    float xmin = nodes_xmin[j];
    float ymin = nodes_ymin[j];
    float width_half = nodes_width[j]/2;
    create_node(xmin, ymin+width_half, width_half);
    nodes_nw[j] = n_nodes-1;
    create_node(xmin+width_half, ymin+width_half, width_half);
    nodes_ne[j] = n_nodes-1;
    create_node(xmin, ymin, width_half);
    nodes_sw[j] = n_nodes-1;
    create_node(xmin+width_half, ymin, width_half);
    nodes_se[j] = n_nodes-1;
}

void eliminate_empty(int j) {
    eliminate_empty_recursive(j, nodes_nw);
    eliminate_empty_recursive(j, nodes_ne);
    eliminate_empty_recursive(j, nodes_sw);
    eliminate_empty_recursive(j, nodes_se);
}

void eliminate_empty_recursive(int j, int children[]) {
    int idx_child = children[j];
    if (idx_child > -1) {
        if (nodes_type[idx_child] == NODE_EXTERNAL) {
            if (nodes_body[idx_child] == -1) {
                children[j] = -1;
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
          ymax = bodies_y[0],
          dx, dy;

    for (i = 1; i < n_bodies; i++) {
        if (bodies_x[i] < xmin) xmin = bodies_x[i];
        if (bodies_x[i] > xmax) xmax = bodies_x[i];
        if (bodies_y[i] < ymin) ymin = bodies_y[i];
        if (bodies_y[i] > ymax) ymax = bodies_y[i];
    }
    dx = xmax-xmin;
    dy = ymax-ymin;

    malloc_nodes();
    create_node(xmin, ymin, dx > dy ? dx : dy);
    for (i = 0; i < n_bodies; i++) {
        insert(0, i);
    }
    eliminate_empty(0);
}

void print_tree(int idx_node, int tab, char *name) {
    int i, idx_body;
    if (idx_node == 0) printf("[TREE]\n");
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

// main loop

void step() {
    int i, epoch = 0;
    float fx, fy;
    clock_t timer_start, timer_end;
    double time_elapsed = 0;
    float time_per_frame;

    while (steps > 0 && epoch < steps) {
        printf("[STEP %d]\n", epoch+1);

        timer_start = clock();

        make_tree();
        for (i = 0; i < n_bodies; i++) {
            compute_force(i, 0, &fx, &fy);
            update_pos(i, fx, fy);
        }

        timer_end = clock();
        time_elapsed += timer_end-timer_start;

        if (is_print_tree) print_tree(0, 0, "root");
        if (is_print_bodies) print_bodies("[UPDATED POSITION]");

        free_nodes();
        epoch += 1;
    }

    time_per_frame = (time_elapsed/(double)CLOCKS_PER_SEC)/epoch;

    printf("\n[TIME]\n");
    printf("Time per iteration: %f s\n", time_per_frame);
    printf("Iteration per second: %f\n", 1/time_per_frame);
}

void handle_sigint(int _) {
    steps = 0;
}

// program

int read_arg(int argc, char *argv[]) {
    int i;
    for (i = 1; i < argc-1; i += 2) {
        if (strcmp(argv[i], "-n_bodies") == 0) {
            if (strcmp(argv[i+1], "data.txt") == 0) {
                if (read_data()) {
                    if (n_bodies < 1) {
                        printf("ERROR: error data.txt\n");
                        return false;
                    }
                } else {
                    return false;
                }
            } else {
                n_bodies = atoi(argv[i+1]);
                if (n_bodies < 1) {
                    printf("ERROR: n_bodies bukan integer atau n_bodies < 1\n");
                    return false;
                }
                malloc_bodies();
                create_body_random();
            }
        } else if (strcmp(argv[i], "-steps") == 0) {
            steps = atoi(argv[i+1]);
            if (steps < 1) {
                printf("ERROR: steps < 1\n");
                return false;
            }
        } else if (strcmp(argv[i], "-print_bodies") == 0) {
            is_print_bodies = true;
        } else if (strcmp(argv[i], "-print_tree") == 0) {
            is_print_tree = true;
        }
    }
    if (n_bodies == 0) {
        printf("USAGE:\n");
        printf("    -n_bodies <data.txt atau integer>    (required)\n");
        printf("        menggunakan body dari data.txt atau\n");
        printf("        generate random sebanyak n_bodies\n");
        printf("    -steps <integer>\n");
        printf("        hanya melakukan iterasi sebanyak steps\n");
        printf("        default: 10\n");
        printf("    -print_bodies 1\n");
        printf("        untuk menampilkan body ke stdout\n");
        printf("    -print_tree 1\n");
        printf("        untuk menampilkan Barnes-Hut tree ke stdout\n");
        return false;
    }
    return true;
}

int read_data() {
    FILE *fp;
    float x, y, mass;
    int i;

    fp = fopen("data.txt", "r");
    if (fp == NULL) {
        printf("Can't open data.txt\n");
        return false;
    }

    n_bodies = 0;
    while (!feof(fp)) {
        if (fscanf(fp, "%f %f %f", &x, &y, &mass) != 3) {
            break;
        }
        n_bodies += 1;
    }
    malloc_bodies();

    i = 0;
    rewind(fp);
    while (!feof(fp)) {
        if (fscanf(fp, "%f %f %f", &x, &y, &mass) != 3) {
            break;
        }
        create_body(i, x, y, mass);
        i += 1;
    }

    fclose(fp);
    return true;
}

// malloc/free

void malloc_bodies() {
    long long int s = n_bodies*sizeof(float);
    bodies_mass = malloc(s);
    bodies_x = malloc(s);
    bodies_y = malloc(s);
    bodies_vx = malloc(s);
    bodies_vy = malloc(s);
}

void malloc_nodes() {
    long long int sf, si;

    n_nodes = n_bodies*16 + 1;
    sf = n_nodes*sizeof(float);
    si = n_nodes*sizeof(int);

    nodes_type = malloc(si);
    nodes_body = malloc(si);
    nodes_mass = malloc(sf);
    nodes_xmass = malloc(sf);
    nodes_ymass = malloc(sf);
    nodes_xmin = malloc(sf);
    nodes_ymin = malloc(sf);
    nodes_width = malloc(sf);
    nodes_nw = malloc(si);
    nodes_ne = malloc(si);
    nodes_sw = malloc(si);
    nodes_se = malloc(si);

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
