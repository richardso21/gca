uniform float iTime;
uniform vec2 iResolution;

//======================================================= Copy-Paste Area Begin =========================================================
// Number of Gaussians
const int NUM_GAUSSIANS = 100;
// Dimensions [x_min, x_max, y_min, y_max]
float dim[4] = float[4](-5.0,5.0,-3.9625000000000004,3.9625000000000004);
// Centers (x, y coordinates)
vec2 gauss_centers[NUM_GAUSSIANS] = vec2[NUM_GAUSSIANS](vec2(5.24, -1.81),vec2(1.19, 3.44),vec2(1.43, -0.00),vec2(1.73, -0.23),vec2(1.78, -3.37),vec2(-0.16, 0.13),vec2(2.03, 2.00),vec2(0.67, -2.75),vec2(4.63, 1.59),vec2(-1.56, -4.15),vec2(4.80, 3.32),vec2(2.67, -2.15),vec2(5.19, 0.42),vec2(3.92, -3.59),vec2(-4.44, 0.32),vec2(-3.60, 0.06),vec2(2.59, 1.45),vec2(0.10, 1.71),vec2(0.74, -1.77),vec2(-1.71, 2.69),vec2(2.13, -2.83),vec2(3.02, 0.67),vec2(0.32, -4.09),vec2(0.07, -2.74),vec2(-2.66, 2.58),vec2(-1.12, 3.77),vec2(1.56, 2.00),vec2(4.04, -0.19),vec2(-0.77, 1.29),vec2(-2.27, 0.19),vec2(3.02, -0.45),vec2(3.87, 1.29),vec2(4.30, 2.09),vec2(-5.07, 0.24),vec2(3.39, -4.14),vec2(3.23, 3.79),vec2(0.85, 2.15),vec2(2.01, 0.67),vec2(-0.14, -1.02),vec2(0.99, -2.18),vec2(-1.45, -0.13),vec2(-3.95, 3.65),vec2(1.33, -3.60),vec2(-4.34, -1.75),vec2(-3.67, 0.16),vec2(-0.30, 1.53),vec2(-1.82, 1.82),vec2(2.98, 0.68),vec2(-1.47, 1.11),vec2(-0.65, 2.52),vec2(1.05, -3.40),vec2(2.74, -1.65),vec2(-0.05, -3.46),vec2(-5.27, 2.58),vec2(-0.97, 2.17),vec2(-4.03, -2.67),vec2(-4.13, 1.52),vec2(0.80, 0.10),vec2(-3.63, -3.96),vec2(-2.33, -1.44),vec2(0.54, 2.58),vec2(-4.95, -2.86),vec2(-3.96, 3.38),vec2(-1.14, -3.23),vec2(1.22, -0.73),vec2(-0.10, 0.91),vec2(-3.43, 2.49),vec2(-2.95, 1.82),vec2(4.08, 0.42),vec2(-1.27, -1.42),vec2(-0.18, -2.00),vec2(-1.57, 1.12),vec2(3.44, 2.26),vec2(3.09, -1.36),vec2(-3.97, -3.05),vec2(-2.70, -1.32),vec2(3.04, 1.47),vec2(2.77, 2.50),vec2(4.00, -2.20),vec2(1.86, -1.19),vec2(5.33, 2.00),vec2(-0.32, 2.49),vec2(-2.61, 1.22),vec2(-5.05, -1.13),vec2(1.56, -2.29),vec2(-2.73, 3.73),vec2(0.69, 1.13),vec2(3.77, -1.16),vec2(4.13, -3.45),vec2(2.07, 2.15),vec2(1.38, 3.20),vec2(3.90, 3.05),vec2(2.52, -2.79),vec2(1.12, 3.27),vec2(-3.99, -1.21),vec2(0.26, 2.33),vec2(4.33, -1.26),vec2(2.31, -3.81),vec2(-1.00, -2.90),vec2(0.35, -3.24));
// Sigmas (scales)
vec2 gauss_sigmas[NUM_GAUSSIANS] = vec2[NUM_GAUSSIANS](vec2(0.27, 0.96),vec2(0.98, 0.58),vec2(-0.00, -0.08),vec2(0.39, 0.61),vec2(-0.00, 0.05),vec2(0.49, 0.53),vec2(-0.00, 0.09),vec2(0.26, 0.34),vec2(0.30, 0.50),vec2(0.42, 0.12),vec2(0.75, 0.66),vec2(0.39, 0.09),vec2(0.46, 0.65),vec2(-0.08, 0.03),vec2(0.54, 0.75),vec2(0.55, 0.16),vec2(-0.00, 0.08),vec2(0.45, 0.36),vec2(0.22, 0.55),vec2(0.94, 0.57),vec2(0.26, 0.04),vec2(0.56, 0.59),vec2(0.52, 0.40),vec2(0.43, 0.19),vec2(0.24, 0.23),vec2(1.01, 0.44),vec2(0.54, 0.60),vec2(0.15, 0.62),vec2(0.36, 0.30),vec2(0.28, 0.73),vec2(0.73, 0.33),vec2(0.35, 0.45),vec2(0.38, 0.51),vec2(0.10, 0.96),vec2(0.18, 0.07),vec2(0.74, 0.60),vec2(0.29, 0.53),vec2(0.52, 0.76),vec2(0.79, 0.40),vec2(0.13, 0.27),vec2(0.53, 0.68),vec2(0.36, 0.47),vec2(0.47, 0.36),vec2(0.84, 0.31),vec2(0.12, 0.85),vec2(0.39, 0.34),vec2(0.39, 0.57),vec2(0.58, 0.68),vec2(0.90, 0.25),vec2(0.07, 0.06),vec2(0.03, 0.06),vec2(0.76, 0.20),vec2(0.00, 0.05),vec2(0.36, 0.88),vec2(0.23, 0.54),vec2(0.28, 0.23),vec2(0.84, 0.62),vec2(0.49, 0.52),vec2(0.72, 0.96),vec2(0.01, -0.05),vec2(0.33, 0.65),vec2(0.55, 1.00),vec2(0.77, 0.92),vec2(1.28, 0.77),vec2(0.57, 0.21),vec2(0.66, 0.29),vec2(0.19, 0.24),vec2(0.55, 0.19),vec2(0.57, 0.34),vec2(0.32, 0.69),vec2(0.24, 0.61),vec2(0.96, 0.43),vec2(0.31, 0.57),vec2(0.30, 0.09),vec2(0.16, 0.04),vec2(0.54, 0.94),vec2(0.42, 0.65),vec2(0.45, 0.54),vec2(0.72, 0.34),vec2(0.80, 0.59),vec2(0.33, 0.88),vec2(0.27, 0.31),vec2(0.14, 0.39),vec2(0.48, 0.63),vec2(0.39, 0.41),vec2(0.14, 0.18),vec2(0.46, 0.39),vec2(0.29, 0.09),vec2(1.00, 0.76),vec2(0.28, 0.38),vec2(0.00, -0.04),vec2(0.49, 0.41),vec2(0.22, 0.55),vec2(0.31, 0.23),vec2(0.32, 0.50),vec2(0.28, 0.31),vec2(0.59, 0.37),vec2(0.40, 0.54),vec2(-0.00, 0.10),vec2(0.13, 0.08));
// Rotation angles (thetas)
float gauss_thetas[NUM_GAUSSIANS] = float[NUM_GAUSSIANS](-0.07,0.11,-0.06,0.65,0.07,-0.47,-0.05,0.56,-0.25,-0.06,-1.02,-0.08,0.42,-0.08,0.35,-0.25,-0.46,0.02,-0.94,-0.10,0.15,-0.01,-0.72,1.14,0.34,0.12,0.48,-0.96,0.55,0.07,0.49,-0.07,-0.22,0.02,0.84,-0.27,-0.06,0.58,0.02,-1.52,-0.04,0.25,-0.41,-0.51,-0.01,0.27,-0.72,0.35,-0.59,-0.31,-0.30,0.06,0.13,0.09,-0.05,-0.31,-1.09,-0.60,-1.06,-0.04,0.66,0.67,0.05,-0.64,0.04,0.42,0.05,0.37,0.25,0.56,-0.64,-0.82,0.54,0.18,-0.18,0.12,0.93,0.21,0.14,-0.18,0.06,-0.40,-0.01,-0.65,0.20,0.39,0.02,0.49,-0.37,0.77,-0.02,0.08,0.88,0.13,0.35,0.11,0.17,-0.98,0.40,-0.02);
// Colors (RGB)
vec3 gauss_colors[NUM_GAUSSIANS] = vec3[NUM_GAUSSIANS](vec3(0.40, 0.42, 0.44),vec3(0.20, 0.26, 0.57),vec3(0.69, 0.30, -0.17),vec3(0.26, 0.34, 0.33),vec3(0.67, 0.05, 0.30),vec3(0.25, 0.33, 0.45),vec3(0.66, 0.29, 0.16),vec3(0.21, 0.24, 0.30),vec3(-0.06, 0.13, 0.42),vec3(0.48, 0.40, 0.46),vec3(0.28, 0.39, 0.50),vec3(0.27, 0.37, 0.36),vec3(0.56, 0.55, 0.50),vec3(0.28, 0.33, 0.15),vec3(0.18, 0.31, 0.50),vec3(0.17, 0.14, -0.01),vec3(0.69, 0.56, 0.44),vec3(0.12, 0.16, 0.33),vec3(0.17, 0.23, 0.38),vec3(0.08, 0.15, 0.44),vec3(0.19, 0.24, 0.26),vec3(0.19, 0.06, 0.17),vec3(0.16, 0.21, 0.24),vec3(0.12, 0.16, 0.16),vec3(0.70, 0.59, -0.14),vec3(0.36, 0.43, 0.49),vec3(0.10, 0.19, 0.38),vec3(0.16, 0.26, 0.39),vec3(0.23, 0.23, 0.29),vec3(0.09, 0.21, 0.41),vec3(0.50, 0.55, 0.46),vec3(-0.01, 0.13, 0.39),vec3(0.64, 0.54, 0.17),vec3(0.82, 0.57, 0.18),vec3(0.26, 0.21, 0.22),vec3(0.21, 0.30, 0.45),vec3(0.25, 0.29, 0.20),vec3(0.29, 0.34, 0.46),vec3(0.39, 0.48, 0.50),vec3(0.18, 0.27, 0.29),vec3(0.73, 0.72, 0.54),vec3(0.55, 0.56, 0.00),vec3(0.20, 0.20, 0.18),vec3(0.07, 0.12, 0.26),vec3(0.27, 0.31, 0.36),vec3(0.19, 0.26, 0.26),vec3(0.51, 0.51, 0.24),vec3(-0.01, 0.22, 0.21),vec3(0.47, 0.36, 0.01),vec3(0.44, 0.20, 0.12),vec3(0.61, 0.29, -0.04),vec3(0.26, 0.33, 0.42),vec3(0.61, 0.15, 0.21),vec3(0.40, 0.49, 0.45),vec3(0.24, 0.30, 0.24),vec3(0.12, 0.15, 0.22),vec3(0.39, 0.46, 0.46),vec3(0.17, 0.27, 0.42),vec3(0.15, 0.17, 0.13),vec3(-0.22, -0.15, 0.11),vec3(-0.04, -0.05, 0.09),vec3(0.20, 0.22, 0.24),vec3(0.01, 0.07, 0.52),vec3(0.11, 0.12, 0.10),vec3(0.41, 0.43, 0.25),vec3(0.10, 0.15, 0.21),vec3(0.10, 0.14, 0.12),vec3(0.20, 0.25, 0.32),vec3(0.56, 0.55, 0.41),vec3(0.11, 0.18, 0.34),vec3(0.15, 0.21, 0.34),vec3(-0.36, -0.16, 0.36),vec3(0.56, 0.54, 0.20),vec3(0.09, 0.14, 0.33),vec3(0.20, 0.21, 0.30),vec3(0.08, 0.09, 0.07),vec3(0.16, 0.20, 0.33),vec3(0.11, 0.21, 0.48),vec3(0.26, 0.33, 0.36),vec3(-0.02, -0.00, 0.17),vec3(0.84, 0.74, 0.42),vec3(0.38, 0.42, 0.36),vec3(0.33, 0.38, 0.42),vec3(0.35, 0.41, 0.36),vec3(0.19, 0.27, 0.35),vec3(0.43, 0.44, 0.29),vec3(0.29, 0.32, 0.35),vec3(0.10, 0.16, 0.26),vec3(0.19, 0.23, 0.27),vec3(0.48, 0.44, 0.12),vec3(0.41, 0.55, -0.14),vec3(0.55, 0.47, 0.05),vec3(0.12, 0.15, 0.19),vec3(0.38, 0.36, -0.21),vec3(0.42, 0.51, 0.54),vec3(0.28, 0.34, 0.27),vec3(0.34, 0.42, 0.64),vec3(0.22, 0.26, 0.32),vec3(0.15, 0.50, 0.31),vec3(0.36, 0.35, 0.43));
//======================================================= Copy-Paste Area End =========================================================

/////////////////////////////////////////////////////
//// Here, you are asked to build the inverse covariance matrix, similar to in the 2DGS_A3_solution.ipynb file.
//// You must create the rotation matrix R, the inverse squared sigma matrix D, and the final inverse covariance matrix.
/////////////////////////////////////////////////////

// This function builds the inverse covariance matrix
mat2 buildSigmaInv(float theta, vec2 sigma)
{
    mat2 cov_mat = mat2(0, 0, 0, 0);

    ///////////
    // BEGINNING OF YOUR CODE.
    //////////
    mat2 R = mat2(cos(theta), sin(theta), -sin(theta), cos(theta));
    vec2 ssigma = abs(sigma) + vec2(0.0001); // safe sigma for division
    mat2 D = mat2(1.0 / (ssigma.x * ssigma.x), 0, 0, 1.0 / (ssigma.y * ssigma.y));
    cov_mat = R * D * transpose(R);

    ///////////
    // END OF YOUR CODE.
    //////////
    return cov_mat;
}

/////////////////////////////////////////////////////
//// Here, you are asked to fill in the necessary components for calculating each gaussian's contribution to the current pixel's color.
//// You must calculate the position of the pixel relative to the gaussian's center, calculate the contribution exponent (pos^T * sigma_inv * pos)
//// and finally calculate the Gaussian function value that will control the contribution of this specific gaussian.
/////////////////////////////////////////////////////

void mainImage(out vec4 fragColor, in vec2 fragCoord)
{
    float aspect = (dim[1] - dim[0]) / (dim[3] - dim[2]) * iResolution.y / iResolution.x;
    vec2 uv = fragCoord.xy / iResolution.xy; // scale to [0, 1]
    // scale from [-1, 1] to [x_min, x_max] and [y_min, y_max]
    uv.x = mix(dim[0], dim[1], uv.x);
    uv.y = mix(dim[2], dim[3], uv.y);
    if(aspect > 1.0)
    {
        uv.y *= aspect;
    }
    else
    {
        uv.x /= aspect;
    }

    vec3 color = vec3(0.0);

    // Draw bounding box
    float edge = 0.01;
    if(uv.x > dim[0] - edge && uv.x < dim[0] + edge ||
        uv.x > dim[1] - edge && uv.x < dim[1] + edge ||
        ((uv.y > dim[2] - edge && uv.y < dim[2] + edge || uv.y > dim[3] - edge && uv.y < dim[3] + edge) && uv.x > dim[0] && uv.x < dim[1]))
    {
        color = vec3(1.0, 1.0, 1.0);
    }

    // Animate the Gaussian centers
    uv += smoothstep(0.0, 1.0, cos(iTime)) * cos(iTime + 100.0 * uv.xy) * 3.;

    for(int i = 0; i < NUM_GAUSSIANS; ++i)
    {
        vec2 center = gauss_centers[i];
        vec2 scale = gauss_sigmas[i];
        float theta = gauss_thetas[i];
        vec3 color_rgb = gauss_colors[i];

        // Build inverse covariance matrix
        mat2 sigma_inv = buildSigmaInv(theta, scale);
        float f_x = 0.;

        ///////////
        // BEGINNING OF YOUR CODE.
        //////////
        vec2 pos = uv - center;
        float a = dot(pos, sigma_inv * pos);
        f_x = exp(-0.5 * a);
        ///////////
        // END OF YOUR CODE.
        //////////

        // Add color contribution
        color += f_x * color_rgb;
    }

    fragColor = vec4(color, 1.0);
}

void main()
{
    mainImage(gl_FragColor, gl_FragCoord.xy);
}