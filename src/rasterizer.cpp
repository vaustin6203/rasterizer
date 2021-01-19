#include "rasterizer.h"
#include "CGL/CGL.h"

using namespace std;

namespace CGL {

RasterizerImp::RasterizerImp(PixelSampleMethod psm, LevelSampleMethod lsm,
                                       size_t width, size_t height,
                                       unsigned int sample_rate) {
  this->psm = psm;
  this->lsm = lsm;
  this->width = width;
  this->height = height;
  this->sample_rate = sample_rate;

  supersample_buffer.resize(width * height * sample_rate, Color::White);
}

// Used by rasterize_point and rasterize_line
void RasterizerImp::fill_pixel(size_t x, size_t y, Color c) {
  // TODO: Task 2: You might need to this function to fix points and lines (such as the black rectangle border in test4.svg)
  // NOTE: You are not required to implement proper supersampling for points and lines
  // It is sufficient to use the same color for all supersamples of a pixel for points and lines (not triangles)
  for (int i = 0; i < this->sample_rate; i++) {
      fill_supersample(x, y, i, c);
  }
}

// Optional helper function to add a sample to the supersample_buffer
void RasterizerImp::fill_supersample(size_t x, size_t y, size_t s, Color c) {
  // TODO: Task 2: You may want to implement this function. Hint: our solution uses one line
  this->supersample_buffer[y * this->width * this->sample_rate + x * this->sample_rate + s] = c;
}

// Rasterize a point: simple example to help you start familiarizing
// yourself with the starter code.
//
void RasterizerImp::rasterize_point(float x, float y, Color color) {
  // fill in the nearest pixel
  int sx = (int)floor(x);
  int sy = (int)floor(y);

  // check bounds
  if (sx < 0 || sx >= width) return;
  if (sy < 0 || sy >= height) return;

  fill_pixel(sx, sy, color);
  return;
}

// Rasterize a line.
void RasterizerImp::rasterize_line(float x0, float y0,
  float x1, float y1,
  Color color) {
  if (x0 > x1) {
    swap(x0, x1); swap(y0, y1);
  }

  float pt[] = { x0,y0 };
  float m = (y1 - y0) / (x1 - x0);
  float dpt[] = { 1,m };
  int steep = abs(m) > 1;
  if (steep) {
    dpt[0] = x1 == x0 ? 0 : 1 / abs(m);
    dpt[1] = x1 == x0 ? (y1 - y0) / abs(y1 - y0) : m / abs(m);
  }

  while (floor(pt[0]) <= floor(x1) && abs(pt[1] - y0) <= abs(y1 - y0)) {
    rasterize_point(pt[0], pt[1], color);
    pt[0] += dpt[0]; pt[1] += dpt[1];
  }
}

//Determines if a point lies inside a triangle
bool inside_test( double x_point, double y_point, double x0, double y0,
                  double x1, double y1, double x2, double y2) {

    //Calculate line for each point on triangle given a point
    double plane01 = -1 * (x_point - x0) * (y1 - y0) + (y_point - y0) * (x1 - x0);
    double plane12 = -1 * (x_point - x1) * (y2 - y1) + (y_point - y1) * (x2 - x1);
    double plane20 = -1 * (x_point - x2) * (y0 - y2) + (y_point - y2) * (x0 - x2);
    double plane10 = -1 * (x_point - x1) * (y0 - y1) + (y_point - y1) * (x0 - x1);
    double plane21 = -1 * (x_point - x2) * (y1 - y2) + (y_point - y2) * (x1 - x2);
    double plane02 = -1 * (x_point - x0) * (y2 - y0) + (y_point - y0) * (x2 - x0);

    return (plane01 >= 0 && plane12 >= 0 && plane20 >= 0) || (plane10 >= 0 && plane21 >= 0 && plane02 >= 0);
}

// Rasterize a triangle.
void RasterizerImp::rasterize_triangle(float x0, float y0,
                                           float x1, float y1,
                                           float x2, float y2,
                                           Color color) {
    // TODO: Task 1: Implement basic triangle rasterization here, no supersampling
    //Find bounding box of triangle
    int min_x = min(min(floor(x0), floor(x1)), floor(x2)) - 1;
    int min_y = min(min(floor(y0), floor(y1)), floor(y2)) - 1;
    int max_x = max(max(ceil(x0), ceil(x1)), ceil(x2)) + 1;
    int max_y = max(max(ceil(y0), ceil(y1)), ceil(y2)) + 1;

    if (min_x < 0 || min_y < 0) {
        return;
    } else if (max_x > this->width || max_y > this->height) {
        return;
    }

    /**
    //If pixel is within the triangle, render it
    float x = min_x;
    for (; max_y >= min_y; max_y --) {
        for (; x <= max_x; x++) {
            if (inside_test(x + 0.5, max_y + 0.5, x0, y0, x1, y1, x2, y2)) {
                fill_pixel(x, max_y, color);
            }
        }
        x = min_x;
    }
*/
  // TODO: Task 2: Update to implement super-sampled rasterization
  int sub_pix_dim = sqrt(this->sample_rate);
  double inc0 = 0.5 * sub_pix_dim;
  double inc1 = 1.f / sub_pix_dim;

  //Traverse up y-axis
  for (int y = min_y -1; y <= max_y + 1; y++) {
      //Traverse up x-axis
      for (int x = min_x -1 ; x <= max_x + 1; x++) {
          int sub_samp_count = 0;
          //Traverse across each subpixel point on y-axis
          for (int sub_y = 0; sub_y < sub_pix_dim; sub_y++) {
              //Traverse across each subpixel point on x-axis
              for (int sub_x = 0; sub_x < sub_pix_dim; sub_x++) {
                  //Is subpixel point inside triangle, add to supersample_buffer
                  if (inside_test(x + inc0 + inc1 * sub_x, y + inc0 + inc1 * sub_y, x0, y0, x1, y1, x2, y2)) {
                      fill_supersample(x, y, sub_samp_count, color);
                  }
                  sub_samp_count += 1;
              }
          }
      }
  }

}

//Calculates a barycentric coordinate
float barycentric_coord(float x, float y,
                        float x_a, float y_a,
                        float x_b, float y_b,
                        float x_c, float y_c) {
    float numerator = -1 * (x - x_b) * (y_c - y_b) + (y - y_b) * (x_c - x_b);
    return numerator / (-1 * (x_a - x_b) * (y_c - y_b) + (y_a - y_b) * (x_c - x_b));
}


void RasterizerImp::rasterize_interpolated_color_triangle(float x0, float y0, Color c0,
                                                          float x1, float y1, Color c1,
                                                          float x2, float y2, Color c2)
{
  // TODO: Task 4: Rasterize the triangle, calculating barycentric coordinates and using them to interpolate vertex colors across the triangle
  // Hint: You can reuse code from rasterize_triangle
    //Find bounding box of triangle
    int min_x = min(min(floor(x0), floor(x1)), floor(x2)) - 1;
    int min_y = min(min(floor(y0), floor(y1)), floor(y2)) - 1;
    int max_x = max(max(ceil(x0), ceil(x1)), ceil(x2)) + 1;
    int max_y = max(max(ceil(y0), ceil(y1)), ceil(y2)) + 1;

    //Check bounds
    if (min_x < 0 || min_y < 0) {
        return;
    } else if (max_x > this->width || max_y > this->height) {
        return;
    }

    //Calculate sample rate increments
    int sub_pix_dim = sqrt(this->sample_rate);
    float inc0 = 0.5 * sub_pix_dim;
    float inc1 = 1.f / sub_pix_dim;

    //Traverse up y-axis
    for (int y = min_y - 1; y <= max_y + 1; y++) {
        //Traverse up x-axis
        for (int x = min_x - 1; x <= max_x + 1; x++) {
            int sub_samp_count = 0;
            //Traverse across each subpixel point on y-axis
            for (int sub_y = 0; sub_y < sub_pix_dim; sub_y++) {
                //Traverse across each subpixel point on x-axis
                for (int sub_x = 0; sub_x < sub_pix_dim; sub_x++) {
                    //Is subpixel point inside triangle, add to supersample_buffer
                    float new_x = x + inc0 + inc1 * sub_x;
                    float new_y = y + inc0 + sub_y * inc1;
                    if (inside_test(new_x, new_y, x0, y0, x1, y1, x2, y2)) {
                        float alpha = barycentric_coord(new_x, new_y, x0, y0, x1, y1, x2, y2);
                        float beta = barycentric_coord(new_x, new_y, x1, y1, x2, y2, x0, y0);
                        float gamma = 1 - alpha - beta;
                        Color new_color = alpha * c0 + beta * c1 + gamma * c2;
                        fill_supersample(x, y, sub_samp_count, new_color);
                    }
                    sub_samp_count += 1;
                }
            }
        }
    }

}

void RasterizerImp::rasterize_textured_triangle(float x0, float y0, float u0, float v0,
    float x1, float y1, float u1, float v1,
    float x2, float y2, float u2, float v2,
    Texture& tex)
{
    // TODO: Task 5: Fill in the SampleParams struct and pass it to the tex.sample function.
    int min_x = min(min(floor(x0), floor(x1)), floor(x2)) - 1;
    int min_y = min(min(floor(y0), floor(y1)), floor(y2)) - 1;
    int max_x = max(max(ceil(x0), ceil(x1)), ceil(x2)) + 1;
    int max_y = max(max(ceil(y0), ceil(y1)), ceil(y2)) + 1;

    //Check bounds
    if (min_x < 0 || min_y < 0) {
        return;
    } else if (max_x > this->width || max_y > this->height) {
        return;
    }

    //Calculate sample rate increments
    int sub_pix_dim = sqrt(this->sample_rate);
    float inc0 = 0.5 * sub_pix_dim;
    float inc1 = 1.f / sub_pix_dim;

    //Fill in SampleParams struct
    SampleParams sp;
    sp.psm = this->psm;
    sp.lsm = this->lsm;

    //Traverse up y-axis
    for (int y = min_y - 1; y <= max_y + 1; y++) {
        //Traverse up x-axis
        for (int x = min_x - 1; x <= max_x + 1; x++) {
            int sub_samp_count = 0;
            //Traverse across each subpixel point on y-axis
            for (int sub_y = 0; sub_y < sub_pix_dim; sub_y++) {
                //Traverse across each subpixel point on x-axis
                for (int sub_x = 0; sub_x < sub_pix_dim; sub_x++) {
                    //Is subpixel point inside triangle, add to supersample_buffer
                    float new_x = x + inc0 + inc1 * sub_x;
                    float new_y = y + inc0 + sub_y * inc1;
                    if (inside_test(new_x, new_y, x0, y0, x1, y1, x2, y2)) {
                        //Calculate textile's uv coords
                        float alpha0 = barycentric_coord(new_x, new_y, x0, y0, x1, y1, x2, y2);
                        float beta0 = barycentric_coord(new_x, new_y, x1, y1, x2, y2, x0, y0);
                        float gamma0 = 1 - alpha0 - beta0;
                        sp.p_uv = Vector2D(alpha0 * u0 + beta0 * u1 + gamma0 * u2, alpha0 * v0 + beta0 * v1 + gamma0 * v2);

                        float alpha1 = alpha0;
                        float beta1 = beta0;
                        float alpha2 = alpha0;
                        float beta2 = beta0;
                        if (inside_test(new_x + 1.0f, new_y, x0, y0, x1, y1, x2, y2)) {
                            alpha1 = barycentric_coord(new_x + 1.0f, new_y, x0, y0, x1, y1, x2, y2);
                            beta1 = barycentric_coord(new_x + 1.0f, new_y, x1, y1, x2, y2, x0, y0);
                        }
                        float gamma1 = 1 - alpha1 - beta1;
                        sp.p_dx_uv = Vector2D(alpha1 * u0 + beta1 * u1 + gamma1 * u2, alpha1 * v0 + beta1 * v1 + gamma1 * v2);

                        if (inside_test(new_x, new_y + 1.0f, x0, y0, x1, y1, x2, y2)) {
                            alpha2 = barycentric_coord(new_x, new_y + 1.0f, x0, y0, x1, y1, x2, y2);
                            beta2 = barycentric_coord(new_x, new_y + 1.0f, x1, y1, x2, y2, x0, y0);
                        }
                        float gamma2 = 1 - alpha2 - beta2;
                        sp.p_dy_uv = Vector2D(alpha2 * u0 + beta2 * u1 + gamma2 * u2, alpha2 * v0 + beta2 * v1 + gamma2 * v2);

                        if (this->lsm == L_ZERO){
                            //If psm = P_LINEAR sample texture color using bilinear interpolation
                            if (this->psm == P_LINEAR) {
                                Color c = tex.sample_bilinear(sp.p_uv, 0);
                                fill_supersample(x, y, sub_samp_count, c);
                                //If psm = P_Nearest sample texture color using nearest neighbor
                            } else if (this->psm == P_NEAREST) {
                                Color c = tex.sample_nearest(sp.p_uv, 0);
                                fill_supersample(x, y, sub_samp_count, c);
                            }
                        } else if (this->lsm == L_NEAREST) {
                            float level = tex.get_level(sp);
                            //If psm = P_LINEAR sample texture color using bilinear interpolation
                            if (this->psm == P_LINEAR) {
                                Color c = tex.sample_bilinear(sp.p_uv, level);
                                fill_supersample(x, y, sub_samp_count, c);
                                //If psm = P_Nearest sample texture color using nearest neighbor
                            } else if (this->psm == P_NEAREST) {
                                Color c = tex.sample_nearest(sp.p_uv, level);
                                fill_supersample(x, y, sub_samp_count, c);
                            }
                        } else if (this->lsm == L_LINEAR) {
                            //finding nearest levels
                            float level0 = tex.get_level(sp);
                            float level1 = floor(level0);
                            float level2 = ceil(level0);

                            //computing weights
                            float diff1 = level0 - level1;
                            float diff2 = level2 - level0;
                            float weight0 = 1 - level0 - level1;
                            float weight1 = 1 - weight0;

                            //If psm = P_LINEAR sample texture color using bilinear interpolation
                            if (this->psm == P_LINEAR) {
                                Color c0 = tex.sample_bilinear(sp.p_uv, level1);
                                Color c1 = tex.sample_bilinear(sp.p_uv, level2);
                                Color c2 = weight0 * c0 + weight1 * c1;
                                fill_supersample(x, y, sub_samp_count, c2);
                                //If psm = P_Nearest sample texture color using nearest neighbor
                            } else if (this->psm == P_NEAREST) {
                                Color c0 = tex.sample_nearest(sp.p_uv, level1);
                                Color c1 = tex.sample_nearest(sp.p_uv, level2);
                                Color c2 = weight0 * c0 + weight1 * c1;
                                fill_supersample(x, y, sub_samp_count, c2);
                            }
                        }
                    }
                    sub_samp_count += 1;
                }
            }
        }
    }
    // TODO: Task 6: Set the correct barycentric differentials in the SampleParams struct.
    // Hint: You can reuse code from rasterize_triangle/rasterize_interpolated_color_triangle

}

void RasterizerImp::set_sample_rate(unsigned int rate) {
  // TODO: Task 2: You may want to update this function for supersampling support

  this->sample_rate = rate;
  this->supersample_buffer.resize(this->width * this->height * rate);
  std::fill(supersample_buffer.begin(), supersample_buffer.end(), Color::White);
}


void RasterizerImp::set_framebuffer_target( unsigned char* rgb_framebuffer,
                                                size_t width, size_t height )
{
  // TODO: Task 2: You may want to update this function for supersampling support

  this->width = width;
  this->height = height;
  this->rgb_framebuffer_target = rgb_framebuffer;
  supersample_buffer.resize(width * height * this->sample_rate);
  std::fill(supersample_buffer.begin(), supersample_buffer.end(), Color::White);
  
}


void RasterizerImp::clear_buffers() {
  // TODO: Task 2: You may want to update this function for supersampling support
  // Hint: With supersampling, you have an additional buffer to take care of
  supersample_buffer.clear();
  supersample_buffer.resize(this->height * this->width * this->sample_rate);
  std::fill(supersample_buffer.begin(), supersample_buffer.end(), Color::White);
  std::fill(rgb_framebuffer_target, rgb_framebuffer_target + 3 * width * height, 255);

}
// This function is called at the end of rasterizing all elements of the
// SVG file.  If you use a supersample buffer to rasterize SVG elements
// for antialising, you could use this call to fill the target framebuffer
// pixels from the supersample buffer data.
//
void RasterizerImp::resolve_to_framebuffer() {
  // TODO: Task 2: You will likely want to update this function for supersampling support
  unsigned int rate = sqrt(this->sample_rate);

  //For each downsampled color, add rgb values to framebuffer
  for (int y = 0; y < this->height; y++) {
      for (int x = 0; x < this->width; x++) {
          Color new_color = Color(0, 0, 0);
          int s = 0;
          for (int y_sub = 0; y_sub < rate; y_sub++) {
              for (int x_sub = 0; x_sub < rate; x_sub++) {
                  //For each supersampled color, downsample color by averaging them
                  new_color += this->supersample_buffer[y * width * this->sample_rate + x * this->sample_rate + s];
                  s++;
              }
          }
          new_color *= (1.f / this->sample_rate);
          rgb_framebuffer_target[3 * (y * width + x)] = (unsigned char) (new_color.r * 255);
          rgb_framebuffer_target[3 * (y * width + x) + 1] = (unsigned char) (new_color.g * 255);
          rgb_framebuffer_target[3 * (y * width + x) + 2] = (unsigned char) (new_color.b * 255);
      }
  }
}

Rasterizer::~Rasterizer() { }

}// CGL
