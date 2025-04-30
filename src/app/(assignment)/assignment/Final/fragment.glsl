#define ARR_COPY(dst, src, n) for (int i = 0; i < n; i++) dst[i] = src[i];

precision highp float;
varying vec2 vUv; // UV (screen) coordinates in [0,1]^2

uniform float iTime;
uniform float iTimeDelta;
uniform float iFrame;
uniform vec2 iResolution;
uniform vec4 iMouse;
uniform sampler2D iChannel0;

const bool CUSTOM = true;

float remap01(float inp, float inp_start, float inp_end) {
    return clamp((inp - inp_start) / (inp_end - inp_start), 0.0, 1.0);
}
float dist_sqr(vec2 a, vec2 b) {
    vec2 diff = a - b;
    return dot(diff, diff);
}

const int MAX_GAUSSIANS = 50;
//======================================================= Copy-Paste Area Begin =========================================================
// Bird
// Number of Gaussians
const int NUM_GAUSSIANS_0 = 20;
// Dimensions [x_min, x_max, y_min, y_max]
float dim_0[4] = float[4](-5.0, 5.0, -5.0, 5.0);
// Centers (x, y coordinates)
const vec2 gauss_centers_0[NUM_GAUSSIANS_0] = vec2[NUM_GAUSSIANS_0](vec2(1.23, 3.42), vec2(-0.43, 2.30), vec2(-2.52, 1.02), vec2(1.71, 0.68), vec2(-2.59, 0.37), vec2(3.49, -0.93), vec2(-1.96, -3.18), vec2(0.01, 2.66), vec2(1.00, -4.21), vec2(2.69, -2.27), vec2(-0.48, 4.33), vec2(0.62, -2.40), vec2(-0.48, -1.44), vec2(1.09, -0.97), vec2(2.49, -0.87), vec2(-1.77, -2.78), vec2(1.44, -0.71), vec2(-1.20, -4.21), vec2(-0.31, -3.21), vec2(-1.25, 0.10));
// Sigmas (scales)
const vec2 gauss_sigmas_0[NUM_GAUSSIANS_0] = vec2[NUM_GAUSSIANS_0](vec2(3.51, 0.43), vec2(4.65, 1.86), vec2(2.42, 3.31), vec2(1.50, 0.85), vec2(1.92, 0.69), vec2(0.77, 0.41), vec2(1.66, 0.80), vec2(2.83, 0.91), vec2(1.46, 0.55), vec2(1.60, 0.95), vec2(2.27, 3.23), vec2(1.48, 2.09), vec2(1.98, 0.91), vec2(0.70, 0.51), vec2(3.27, 1.53), vec2(2.82, 2.23), vec2(2.46, 1.65), vec2(1.58, 2.62), vec2(0.63, 1.56), vec2(4.09, 1.14));
// Rotation angles (thetas)
const float gauss_thetas_0[NUM_GAUSSIANS_0] = float[NUM_GAUSSIANS_0](-0.79, -1.42, -0.83, -0.15, 1.19, 0.06, 2.20, -0.73, 0.10, 0.12, 1.71, -2.84, -1.38, -0.43, -0.32, 2.23, 3.14, 0.77, 1.88, 1.18);
// Colors (RGB)
const vec3 gauss_colors_0[NUM_GAUSSIANS_0] = vec3[NUM_GAUSSIANS_0](vec3(0.84, -0.03, 0.19), vec3(-1.46, 0.18, -0.33), vec3(0.33, 0.08, 0.02), vec3(1.33, 0.23, 0.47), vec3(0.94, -0.06, 0.18), vec3(0.87, 0.92, 1.09), vec3(1.22, 0.13, 0.38), vec3(1.11, -0.10, 0.27), vec3(1.35, 1.47, 1.13), vec3(1.79, 1.11, 0.44), vec3(0.14, -0.10, 0.05), vec3(0.31, -1.21, -0.39), vec3(1.63, 0.31, 0.59), vec3(1.63, 1.35, 1.51), vec3(-0.57, -1.11, -0.22), vec3(-0.42, 0.24, 0.31), vec3(-0.29, 0.99, -0.17), vec3(-0.78, -0.63, -0.78), vec3(0.50, 1.30, 0.92), vec3(1.01, -0.01, 0.20));

// Pig
// Number of Gaussians
const int NUM_GAUSSIANS_1 = 30;
// Dimensions [x_min, x_max, y_min, y_max]
const float dim_1[4] = float[4](-5.0, 5.0, -5.0, 5.0);
// Centers (x, y coordinates)
const vec2 gauss_centers_1[NUM_GAUSSIANS_1] = vec2[NUM_GAUSSIANS_1](vec2(-0.83, 0.73), vec2(4.04, 0.55), vec2(-0.64, -1.24), vec2(3.33, -0.76), vec2(1.08, -0.96), vec2(-3.55, 1.41), vec2(-3.79, -0.83), vec2(-0.82, -1.20), vec2(1.70, 2.09), vec2(2.21, 2.42), vec2(-3.00, 3.64), vec2(-0.39, 1.87), vec2(-3.75, -0.83), vec2(3.02, -0.75), vec2(-1.63, -0.33), vec2(-0.40, -2.11), vec2(2.98, 1.36), vec2(-0.86, 3.34), vec2(0.66, 4.30), vec2(-0.95, 0.84), vec2(-2.49, -3.26), vec2(1.65, -3.44), vec2(0.65, 0.50), vec2(-0.59, 2.94), vec2(-3.76, -0.09), vec2(0.68, 0.84), vec2(-2.59, 1.03), vec2(3.04, -0.76), vec2(-0.53, -0.55), vec2(2.29, 2.43));
// Sigmas (scales)
const vec2 gauss_sigmas_1[NUM_GAUSSIANS_1] = vec2[NUM_GAUSSIANS_1](vec2(0.09, -0.63), vec2(0.86, 0.17), vec2(0.41, 0.50), vec2(0.17, 0.17), vec2(-0.59, 0.32), vec2(0.28, 1.44), vec2(0.55, 0.66), vec2(-0.58, 0.77), vec2(1.45, 0.42), vec2(-0.14, -0.61), vec2(0.48, 0.31), vec2(0.65, 1.39), vec2(0.26, 0.31), vec2(0.78, 0.71), vec2(0.48, -0.86), vec2(0.57, 1.29), vec2(0.27, 0.08), vec2(-1.22, 0.13), vec2(0.43, 0.54), vec2(2.58, 0.34), vec2(1.65, 0.50), vec2(1.69, 0.59), vec2(0.43, -0.08), vec2(0.35, 1.07), vec2(3.64, 1.00), vec2(0.01, -0.00), vec2(1.30, 0.73), vec2(1.05, 0.88), vec2(-0.53, 1.40), vec2(0.88, 0.18));
// Rotation angles (thetas)
const float gauss_thetas_1[NUM_GAUSSIANS_1] = float[NUM_GAUSSIANS_1](-1.47, -1.14, 0.35, -0.84, -1.40, -0.61, -0.20, -2.16, -0.63, 1.00, 0.48, 1.46, -0.20, 2.11, -0.77, -1.56, -0.16, 0.18, -1.34, 0.04, -0.40, 0.27, -0.47, -1.48, 0.80, -0.01, 1.13, -0.88, 2.12, -0.58);
// Colors (RGB)
const vec3 gauss_colors_1[NUM_GAUSSIANS_1] = vec3[NUM_GAUSSIANS_1](vec3(-0.23, -0.46, -0.24), vec3(0.31, 0.65, 0.19), vec3(2.01, 2.50, -0.24), vec3(-1.65, -1.38, -1.71), vec3(0.57, 0.67, -0.08), vec3(0.32, 0.67, 0.22), vec3(2.32, 2.29, 2.50), vec3(-2.14, -2.70, 0.23), vec3(0.31, 0.64, 0.20), vec3(-0.67, -1.35, -0.48), vec3(0.50, 1.09, 0.33), vec3(0.52, 1.02, 0.41), vec3(-2.26, -2.30, -2.41), vec3(0.81, -1.23, 1.33), vec3(0.95, 1.21, -0.08), vec3(1.10, 1.52, -0.03), vec3(-0.24, -0.49, -0.23), vec3(0.31, 0.59, 0.24), vec3(0.53, 1.08, 0.39), vec3(0.16, 0.40, 0.19), vec3(0.44, 1.01, 0.34), vec3(0.50, 1.08, 0.36), vec3(-0.33, -0.48, -0.16), vec3(0.40, 0.72, 0.32), vec3(-0.33, -0.39, -0.38), vec3(-0.05, 0.72, 0.04), vec3(0.61, 0.98, 0.53), vec3(0.52, 2.30, 0.06), vec3(0.86, 1.21, 0.01), vec3(0.88, 1.81, 0.62));
//======================================================= Copy-Paste Area End =========================================================

// This function builds the inverse covariance matrix
mat2 buildSigmaInv(float theta, vec2 sigma) {
    mat2 cov_mat = mat2(0, 0, 0, 0);

    mat2 R = mat2(cos(theta), sin(theta), -sin(theta), cos(theta));
    vec2 ssigma = abs(sigma) + vec2(0.0001); // safe sigma for division
    mat2 D = mat2(1.0 / (ssigma.x * ssigma.x), 0, 0, 1.0 / (ssigma.y * ssigma.y));
    cov_mat = R * D * transpose(R);

    return cov_mat;
}

// ------------------------------------------------------------
// Particle structure
struct Particle {
    vec2 pos;
    vec2 pos_prev;
    vec2 vel;
    float inv_mass;
    bool is_fixed;
};

// Simulation constants
const float damp = 0.4;
const float collision_dist = 0.2;
const float ground_collision_dist = 0.1;
const vec2 gravity = vec2(0.0, -1);

// Define n_rope rope particles and add one extra "mouse particle".
const int MAX_PARTICLES = 30;
const int MAX_SPRINGS = 30;
const int STARTING_SPRINGS = 3;
//0: mouse particle
//1...5: rope particles
int n_particles;
Particle particles[MAX_PARTICLES];

int nearest_particle(vec2 p) {
    int idx = 1;
    float min_dist = 1e9;
    for(int i = 1; i < n_particles; i++) {
        float d = dist_sqr(p, particles[i].pos);
        if(d < min_dist) {
            min_dist = d;
            idx = i;
        }
    }
    return idx;
}

// ------------------------------------------------------------
// Spring structure
struct Spring {
    int a;
    int b;
    float restLength;
    float inv_stiffness;
};
// Create springs between adjacent rope particles (n_rope-1 springs)
// and one spring connecting the last rope particle and the mouse particle.
Spring springs[MAX_SPRINGS];
int n_springs;
int selected_particle = -1;
int current_add_particle = -1;

Spring add_spring(int a, int b, float inv_stiffness) {
    Spring s;
    s.a = a;
    s.b = b;
    s.restLength = length(particles[a].pos - particles[b].pos);
    s.inv_stiffness = inv_stiffness;
    return s;
}

const int initial_particles = CUSTOM ? 9 : 6;

void init_state(void) {
    if(!CUSTOM) {
        n_particles = 6;
        n_springs = 5;

        //particle 0 is the mouse particle and will be set later
        particles[1].pos = vec2(-0.6, 0.5);
        particles[1].vel = vec2(0.0);
        particles[2].pos = vec2(-0.3, 0.5);
        particles[2].vel = vec2(0.0);
        particles[3].pos = vec2(-0, 0.5);
        particles[3].vel = vec2(0.0);
        particles[4].pos = vec2(0.3, 0.5);
        particles[4].vel = vec2(0.0);
        particles[5].pos = vec2(0.6, 0.5);
        particles[5].vel = vec2(0.0);

        // Springs between adjacent rope particles
        //spring 0 is the mouse particle to the first rope particle
        springs[1] = add_spring(1, 2, 1.0 / 100.0); // first to second rope particle
        springs[2] = add_spring(2, 3, 1.0 / 100.0); // second to third rope particle
        springs[3] = add_spring(3, 4, 1.0 / 100.0); // third to fourth rope particle
        springs[4] = add_spring(4, 5, 1.0 / 100.0); // fourth to fifth rope particle
    } else {
        // TODO: add initial states
        n_particles = 5;
        
        // Bird
        particles[1].pos = vec2(-1.25, -.1);
        particles[1].vel = vec2(0.);

        // Slingshot
        particles[2].pos = vec2(-1.5, 0.2);
        particles[2].vel = vec2(0.);
        particles[3].pos = vec2(-1., -.4);
        particles[3].vel = vec2(0.);

        // Pig
        particles[4].pos = vec2(1.1, -0.4);

        n_springs = STARTING_SPRINGS;
        springs[1] = add_spring(1, 2, 1.0 / 1000.0);
        springs[2] = add_spring(1, 3, 1.0 / 1000.0);

        // Blocks
        vec2 centers[4] = vec2[](vec2(0.4, -0.3), vec2(0.8, -0.3), vec2(0.6, 0.1), vec2(1.4, -0.3));
        for (int i = 0; i < 4; i++) {
            particles[n_particles].pos = centers[i] + vec2(0.1, 0.1);
            particles[n_particles+1].pos = centers[i] + vec2(-0.1, 0.1);
            particles[n_particles+2].pos = centers[i] + vec2(0.1, -0.1);
            particles[n_particles+3].pos = centers[i] + vec2(-0.1, -0.1);

            springs[n_springs] = add_spring(n_particles, n_particles+1, 1.0 / 1000.0);
            springs[n_springs+1] = add_spring(n_particles+1, n_particles+2, 1.0 / 1000.0);
            springs[n_springs+2] = add_spring(n_particles+2, n_particles+3, 1.0 / 1000.0);
            springs[n_springs+3] = add_spring(n_particles+3, n_particles, 1.0 / 1000.0);
            springs[n_springs+4] = add_spring(n_particles, n_particles+2, 1.0 / 1000.0); 
            springs[n_springs+5] = add_spring(n_particles+1, n_particles+3, 1.0 / 1000.0); 
            
            n_particles = n_particles + 4;
            n_springs = n_springs + 6;
        }
    }
    current_add_particle = initial_particles;
}

vec2 screen_to_xy(vec2 coord) {
    return (coord - 0.5 * iResolution.xy) * 2.0 / iResolution.y;
}

bool is_initializing() {
    return iTime < 0.06 || iFrame < 2.;
}

// Load rope particles from the previous frame and update the mouse particle.
void load_state() {
    //0,0: (num_particles, num_springs, selected_particle)

    vec4 data = texelFetch(iChannel0, ivec2(0, 0), 0);
    n_particles = int(data.x);
    n_springs = int(data.y);
    selected_particle = int(data.z);
    current_add_particle = int(data.w);

    //initialize mouse particle
    {
        int mouse_idx = 0;
        particles[mouse_idx].pos = screen_to_xy(iMouse.xy);
        particles[mouse_idx].vel = vec2(0.0);
        particles[mouse_idx].inv_mass = 0.0; // fixed particle
        particles[mouse_idx].is_fixed = true;
    }
    // Load other particles
    for(int i = 1; i < n_particles; i++) {
        vec4 data = texelFetch(iChannel0, ivec2(i, 0), 0);
        particles[i].pos = data.xy;
        particles[i].vel = data.zw;
        particles[i].inv_mass = 1.0; // all particles have mass 1.0
        particles[i].is_fixed = false;

        if(!CUSTOM && (i == 1 || i == 5) || CUSTOM && (i == 2 || i == 3)) {
            particles[i].inv_mass = 0.0; // fixed particles at the ends of the rope
            particles[i].is_fixed = true; // make sure the first and last particles are fixed
        }
    }

    //select nearest particle to mouse
    if(iMouse.z == 1.) {
        if(selected_particle == -1) {
            selected_particle = nearest_particle(particles[0].pos);
        }
    } else {
        selected_particle = -1;
    }

    if(iMouse.z == 2.) {
        particles[current_add_particle].pos = screen_to_xy(iMouse.xy); // update the position of the selected particle
        particles[current_add_particle].vel = vec2(0.0); // reset velocity to zero when mouse is released
        particles[current_add_particle].inv_mass = 1.0; // make sure the selected particle is fixed
        particles[current_add_particle].is_fixed = false; // make sure the selected particle is fixed
        if(current_add_particle >= n_particles) {
            // If we reach the maximum number of particles, reset to the first available index.
            n_particles = current_add_particle + 1; // skip the mouse particle at index 0
        }
        current_add_particle++;
        if(current_add_particle >= MAX_PARTICLES) {
            current_add_particle = initial_particles;
        }
    }

    //load springs
    springs[0] = Spring(0, selected_particle, 0.0, 1.0 / 1000.0); // mouse particle to first rope particle
    for(int i = 1; i < n_springs; i++) {
        vec4 data = texelFetch(iChannel0, ivec2(i, 1), 0);
        springs[i].a = int(data.x);
        springs[i].b = int(data.y);
        springs[i].restLength = data.z;
        springs[i].inv_stiffness = data.w;
    }
}

float spring_constraint(Spring s) {
    return length(particles[s.a].pos - particles[s.b].pos) - s.restLength;
}

vec2 spring_constraint_gradient(vec2 a, vec2 b) {
    return (length(a - b) == 0.0) ? vec2(0.0) : (a - b) / length(a - b);
}

// Compute the gradient of the spring constraint with respect to a given particle.
vec2 spring_constraint_grad(Spring s, int particle_idx) {
    float sgn = (particle_idx == s.a) ? 1.0 : -1.0;
    return sgn * spring_constraint_gradient(particles[s.a].pos, particles[s.b].pos);
}

void solve_spring(Spring s, float dt) {
    float numer = 0.;
    float denom = 0.;

    numer = -spring_constraint(s);
    vec2 grad_a = spring_constraint_grad(s, s.a);
    vec2 grad_b = spring_constraint_grad(s, s.b);
    denom += particles[s.a].inv_mass * pow(length(grad_a), 2.);
    denom += particles[s.b].inv_mass * pow(length(grad_b), 2.);

    // PBD if you comment out the following line
    denom += s.inv_stiffness / (dt * dt);

    if(denom == 0.0)
        return;
    float lambda = numer / denom;
    particles[s.a].pos += lambda * particles[s.a].inv_mass * grad_a;
    particles[s.b].pos += lambda * particles[s.b].inv_mass * grad_b;
}

float collision_constraint(vec2 a, vec2 b, float collision_dist) {
    float dist = length(a - b);
    if(dist < collision_dist) {
        return dist - collision_dist;
    } else {
        return 0.0;
    }
}

vec2 collision_constraint_gradient(vec2 a, vec2 b, float collision_dist) {
    float dist = length(a - b);
    if(dist < collision_dist) {
        return (a - b) / dist;
    } else {
        return vec2(0.0, 0.0);
    }
}

void solve_collision_constraint(int i, int j, float collision_dist, float dt) {
    float numer = 0.0;
    float denom = 0.0;

    numer = -collision_constraint(particles[i].pos, particles[j].pos, collision_dist);
    vec2 grad = collision_constraint_gradient(particles[i].pos, particles[j].pos, collision_dist);
    denom += particles[i].inv_mass * pow(length(grad), 2.);
    denom += particles[j].inv_mass * pow(length(grad), 2.);

    //PBD if you comment out the following line, which is faster
    denom += (1. / 1000.) / (dt * dt);

    if(denom == 0.0)
        return;
    float lambda = numer / denom;
    particles[i].pos += lambda * particles[i].inv_mass * grad;
    particles[j].pos -= lambda * particles[j].inv_mass * grad;
}

float phi(vec2 p) {
    const float PI = 3.14159265359;
    //let's do sin(x)+0.5
    if(!CUSTOM)
        return p.y - (0.1 * sin(p.x * 2. * PI) - 0.5);
    else
        return p.y + 0.5;
}

float ground_constraint(vec2 p, float ground_collision_dist) {
    if(phi(p) < ground_collision_dist) {
        return phi(p) - ground_collision_dist;
    } else {
        return 0.0;
    }
}

vec2 ground_constraint_gradient(vec2 p, float ground_collision_dist) {
    if(phi(p) < ground_collision_dist) {
        const float PI = 3.14159265359;
        return vec2((-0.1 * 2. * PI * cos(2. * PI * p.x)), 1.0);
    } else {
        return vec2(0.0, 0.0);
    }
}

void solve_ground_constraint(int i, float ground_collision_dist, float dt) {
    float numer = 0.0;
    float denom = 0.0;

    vec2 grad = ground_constraint_gradient(particles[i].pos, ground_collision_dist);
    numer = -ground_constraint(particles[i].pos, ground_collision_dist);
    denom += particles[i].inv_mass * pow(length(grad), 2.);

    //PBD if you comment out the following line, which is faster
    denom += (1. / 1000.) / (dt * dt);

    if(denom == 0.0)
        return;
    float lambda = numer / denom;

    // Add friction to ground
    particles[i].pos += lambda * particles[i].inv_mass * vec2(grad.x * .1, grad.y);
}

float boundary_constraint(vec2 p, vec2 bound) {
    float dist_x = max(abs(p.x) - bound.x, 0.0);
    float dist_y = max(abs(p.y) - bound.y, 0.0);
    return max(dist_x, dist_y);
}
vec2 boundary_constraint_gradient(vec2 p, vec2 bound) {
    vec2 grad = vec2(0.0);
    if(abs(p.x) > bound.x) {
        grad.x = sign(p.x);
    }
    if(abs(p.y) > bound.y) {
        grad.y = sign(p.y);
    }
    return grad;
}
void solve_boundary_constraint(int i, vec2 bound, float dt) {
    float numer = -boundary_constraint(particles[i].pos, bound);
    vec2 grad = boundary_constraint_gradient(particles[i].pos, bound);
    float denom = particles[i].inv_mass * dot(grad, grad);
    denom += (1.0 / 1000.0) / (dt * dt);

    if(denom == 0.0)
        return;
    float lambda = numer / denom;
    particles[i].pos += lambda * particles[i].inv_mass * grad;
}

void solve_constraints(float dt) {
    if(iMouse.z == 1.) {
        solve_spring(springs[0], dt); // mouse particle to first rope particle
    }

    // Solve all constraints
    for(int i = 1; i < n_springs; i++) {
        solve_spring(springs[i], dt);
    }
    for(int i = 1; i < n_particles; i++) {
        for(int j = i + 1; j < n_particles; j++) {
            solve_collision_constraint(i, j, collision_dist, dt);
        }
    }
    for(int i = 1; i < n_particles; i++) {
        solve_ground_constraint(i, ground_collision_dist, dt);
        if(CUSTOM)
            solve_boundary_constraint(i, vec2(1.5, 0.9), dt);
    }
}

float dist_to_segment(vec2 p, vec2 a, vec2 b) {
    vec2 pa = p - a;
    vec2 ba = b - a;
    // Compute the projection factor and clamp it between 0 and 1.
    float h = clamp(dot(pa, ba) / dot(ba, ba), 0.0, 1.0);
    // Return the distance from p to the closest point on the segment.
    return length(pa - h * ba);
}

vec3 gaussianImage(in vec2 fragCoord, in int image_idx) {
    // select correct gaussian parameters to draw
    int NUM_GAUSSIANS;
    float dim[4];
    vec2 gauss_centers[MAX_GAUSSIANS];
    vec2 gauss_sigmas[MAX_GAUSSIANS];
    float gauss_thetas[MAX_GAUSSIANS];
    vec3 gauss_colors[MAX_GAUSSIANS];

    if(image_idx == 0) {
        NUM_GAUSSIANS = NUM_GAUSSIANS_0;
        dim = dim_0;
        ARR_COPY(gauss_centers, gauss_centers_0, NUM_GAUSSIANS);
        ARR_COPY(gauss_sigmas, gauss_sigmas_0, NUM_GAUSSIANS);
        ARR_COPY(gauss_thetas, gauss_thetas_0, NUM_GAUSSIANS);
        ARR_COPY(gauss_colors, gauss_colors_0, NUM_GAUSSIANS);
    } else if(image_idx == 1) {
        NUM_GAUSSIANS = NUM_GAUSSIANS_1;
        dim = dim_1;
        ARR_COPY(gauss_centers, gauss_centers_1, NUM_GAUSSIANS);
        ARR_COPY(gauss_sigmas, gauss_sigmas_1, NUM_GAUSSIANS);
        ARR_COPY(gauss_thetas, gauss_thetas_1, NUM_GAUSSIANS);
        ARR_COPY(gauss_colors, gauss_colors_1, NUM_GAUSSIANS);
    }

    // change scaling of the gaussian according to the screen's resolution
    float scale = 2.5 / iResolution.y;
    dim[0] = dim[0] * (1. / scale);  // x_min
    dim[1] = dim[1] * (1. / scale);  // x_max
    dim[2] = dim[2] * (1. / scale);  // y_min
    dim[3] = dim[3] * (1. / scale);  // y_max

    float aspect = (dim[1] - dim[0]) / (dim[3] - dim[2]) * iResolution.y / iResolution.x;
    vec2 uv = (fragCoord.xy) / (iResolution.xy); // scale to [0, 1]

    // scale from [-1, 1] to [x_min, x_max] and [y_min, y_max]
    uv.x = mix(dim[0], dim[1], uv.x);
    uv.y = mix(dim[2], dim[3], uv.y);

    if(aspect > 1.0) {
        uv.y *= aspect;
    } else {
        uv.x /= aspect;
    }

    vec3 color = vec3(0.0);

    for(int i = 0; i < NUM_GAUSSIANS; ++i) {
        vec2 center = gauss_centers[i];
        vec2 scale = gauss_sigmas[i];
        float theta = gauss_thetas[i];
        vec3 color_rgb = gauss_colors[i];

        // Build inverse covariance matrix
        mat2 sigma_inv = buildSigmaInv(theta, scale);
        float f_x = 0.;

        vec2 pos = uv - center;
        float a = dot(pos, sigma_inv * pos);
        f_x = exp(-0.5 * a);

        // Add color contribution
        color += f_x * color_rgb;
    }

    return color;
}

vec3 render_scene(vec2 pixel_xy) {
    float phi = phi(pixel_xy);
    vec3 col;
    if(phi < 0.0) {
        col = vec3(122, 183, 0) / 255.; // ground color
    } else {
        // col = vec3(229, 242, 250) / 255.; // background color
        col = vec3(0.);
    }

    float pixel_size = 2.0 / iResolution.y;

    // If still initializing, return the background color.
    if(is_initializing()) {
        return col;
    }

    // Render rope particles
    {
        float min_dist = 1e9;
        int min_dist_index = 0;

        if(iMouse.z == 1.) {
            min_dist = dist_sqr(pixel_xy, particles[0].pos);
        }

        for(int i = 1; i < n_particles; i++) {
            float dist = dist_sqr(pixel_xy, particles[i].pos);
            if(dist < min_dist) {
                min_dist = dist;
                min_dist_index = i;
            }
        }

        min_dist = sqrt(min_dist);

        const float radius = 0.1;
        if(min_dist_index == 1) { // render bird
            vec2 particle_pos = particles[1].pos;
            // Transform pixel coordinates relative to particle position
            vec2 gaussian_coord = (pixel_xy - particle_pos) * 20.0 + iResolution.xy / 2.0;

            // Only render within a certain radius of the particle
            float dist_to_particle = length(pixel_xy - particle_pos);
            if(dist_to_particle < 0.1) { // Adjust radius as needed
                vec3 gaussian_color = gaussianImage(gaussian_coord, 0);
                // Blend gaussian color with background
                col = mix(col, gaussian_color, 1.0);
            }
        } else if(min_dist_index == 4) { // render bird
            vec2 particle_pos = particles[4].pos;
            // Transform pixel coordinates relative to particle position
            vec2 gaussian_coord = (pixel_xy - particle_pos) * 20.0 + iResolution.xy / 2.0;

            // Only render within a certain radius of the particle
            float dist_to_particle = length(pixel_xy - particle_pos);
            if(dist_to_particle < 0.1) { // Adjust radius as needed
                vec3 gaussian_color = gaussianImage(gaussian_coord, 1);
                // Blend gaussian color with background
                col = mix(col, gaussian_color, 1.0);
            }
        // Only render regular particles for non-gaussian particles
        } 
        else
            col = mix(col, vec3(180, 164, 105) / 255., remap01(min_dist, radius, radius - pixel_size));
    }

    // Render All springs
    {
        float min_dist = 1e9;

        if(iMouse.z == 1.) {
            min_dist = dist_to_segment(pixel_xy, particles[0].pos, particles[selected_particle].pos);
        }

        vec3 color = vec3(14, 105, 146) / 255.;
        float opacity = 0.25;
        float thickness = 0.01;
        for(int i = 1; i < STARTING_SPRINGS; i++) {
            int a = springs[i].a;
            int b = springs[i].b;
            min_dist = min(min_dist, dist_to_segment(pixel_xy, particles[a].pos, particles[b].pos));
        }

        for (int i = STARTING_SPRINGS; i < n_springs; i++) {
            int a = springs[i].a;
            int b = springs[i].b;
            if (dist_to_segment(pixel_xy, particles[a].pos, particles[b].pos) < min_dist){
                min_dist = dist_to_segment(pixel_xy, particles[a].pos, particles[b].pos);
                color = vec3(180, 164, 105) / 255.;
                thickness = .1;
                opacity = 1.;
            }
        }

        col = mix(col, color, opacity * remap01(min_dist, thickness, thickness - pixel_size));
    }

    // col.z = 1.0;
    return col;
}

vec4 output_color(vec2 pixel_ij) {
    int i = int(pixel_ij.x);
    int j = int(pixel_ij.y);

    if(j == 0) {
        // (0,0): (num_particles, num_springs, selected_particle)
        if(i == 0) {
            return vec4(float(n_particles), float(n_springs), float(selected_particle), float(current_add_particle));
        } else if(i < n_particles) {
            //a particle
            return vec4(particles[i].pos, particles[i].vel);
        } else {
            return vec4(0.0, 0.0, 0.0, 1.0);
        }
    } else if(j == 1) {
        if(i < n_springs) {
            return vec4(float(springs[i].a), float(springs[i].b), springs[i].restLength, springs[i].inv_stiffness);
        } else {
            return vec4(0.0, 0.0, 0.0, 1.0);
        }
    } else {
        vec2 pixel_xy = screen_to_xy(pixel_ij);
        vec3 color = render_scene(pixel_xy);
        return vec4(color, 1.0);
    }
}

// ------------------------------------------------------------
// Main function
void main() {
    vec2 pixel_ij = vUv * iResolution.xy;
    int pixel_i = int(pixel_ij.x);
    int pixel_j = int(pixel_ij.y);

    if(is_initializing()) {
        init_state();
    } else {
        load_state();
        if(pixel_j == 0) {
            if(pixel_i >= n_particles)
                return;

            float actual_dt = min(iTimeDelta, 0.02);
            const int n_steps = 5;
            float dt = actual_dt / float(n_steps);

            for(int i = 0; i < n_steps; i++) {
                // Update rope particles only; skip updating the mouse particle since it's fixed.
                for(int j = 0; j < n_particles; j++) {
                    if(!particles[j].is_fixed)
                        particles[j].vel += dt * gravity;
                    particles[j].vel *= exp(-damp * dt);
                    particles[j].pos_prev = particles[j].pos;
                    particles[j].pos += dt * particles[j].vel;
                }
                solve_constraints(dt);
                // Update velocities for rope particles only.
                for(int j = 0; j < n_particles; j++) {
                    if(!particles[j].is_fixed) {
                        particles[j].vel = (particles[j].pos - particles[j].pos_prev) / dt;
                    }
                }
                // Keep the mouse particle fixed by reassigning its position each step.
                int mouse_idx = 0;
                particles[mouse_idx].pos = screen_to_xy(iMouse.xy);
            }
        }
    }

    gl_FragColor = output_color(pixel_ij);
}
