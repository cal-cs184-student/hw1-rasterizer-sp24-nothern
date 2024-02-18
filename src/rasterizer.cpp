#include "rasterizer.h"
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

    sample_buffer.resize(width * height * sample_rate, Color::White);
  }

  // Used by rasterize_point and rasterize_line
  void RasterizerImp::fill_pixel(size_t x, size_t y, Color c) {
    // TODO: Task 2: You might need to this function to fix points and lines (such as the black rectangle border in test4.svg)
    // NOTE: You are not required to implement proper supersampling for points and lines
    // It is sufficient to use the same color for all supersamples of a pixel for points and lines (not triangles)

    double rate_div = sqrt(sample_rate);
    sample_buffer[y * width * rate_div + x] = c;
  }

  // Rasterize a point: simple example to help you start familiarizing
  // yourself with the starter code.
  //
  void RasterizerImp::rasterize_point(float x, float y, Color color) {
    // fill in the nearest pixel
    int sx = (int)floor(x);
    int sy = (int)floor(y);

    // adjust
    double rate_div = sqrt(sample_rate);

    

    // check bounds
    if (sx < 0 || sx >= rate_div * width) return;
    if (sy < 0 || sy >= rate_div * height) return;

    /*float biass = 0.5f / rate_div;*/

    /*float sample_step = 1.0f / rate_div*/;

    float left_border_x = rate_div * sx;
    float left_border_y = rate_div * sy;

    float right_border_x = rate_div * (sx + 1);
    float right_border_y = rate_div * (sy + 1);

    for (float temp_x = left_border_x; temp_x < right_border_x; temp_x++) {
        for (float temp_y = left_border_y; temp_y < right_border_y; temp_y++) {
            //float sample_x = temp_x + 0.5f;
            //float sample_y = temp_y + 0.5f;

            bool inside = true;

            if (inside) {
                fill_pixel(static_cast<size_t>(temp_x), static_cast<size_t>(temp_y), color);
            }
        }
    }

    /*fill_pixel(sx, sy, color);*/
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

  // Rasterize a triangle.
  void RasterizerImp::rasterize_triangle(float x0, float y0,
    float x1, float y1,
    float x2, float y2,
    Color color) {
    // TODO: Task 1: Implement basic triangle rasterization here, no supersampling

    // the general process is:
    // 1 find the 4 edges for the whole rectangle to cover the triangle
    // 2 3-line-test from lecture 2
    // 3 iteration plot
    // 4 optimize
    
    // refer to Ed#81bea for test6
    // If area is negative, the vertices are in clockwise order, so swap them
    float direction = (x1 - x0) * (y2 - y0) - (x2 - x0) * (y1 - y0);
    if(direction < 0) {
        swap(x1, x2);
        swap(y1, y2);
    }

    // find many white line inside the image. refer to ed81a, use floor and ceil to solve
    float leftx = floor(min({x0, x1, x2}));
    float lowy = floor(min({y0, y1, y2}));

    float rightx = ceil(max({x0, x1, x2}));
    float highy = ceil(max({y0, y1, y2}));
    
    /*for(float x = leftx; x <= rightx; ++x) {
        for(float y = lowy; y <= highy; ++y) {
            float sample_x = x + 0.5f;
            float sample_y = y + 0.5f;

            bool inside = ((x1 - x0) * (sample_y - y0) - (y1 - y0) * (sample_x - x0)) >= 0 && ((x2 - x1) * (sample_y - y1) - (y2 - y1) * (sample_x - x1)) >= 0 && ((x0 - x2) * (sample_y - y2) - (y0 - y2) * (sample_x - x2)) >= 0;

            if(inside) {
                fill_pixel(static_cast<size_t>(sample_x), static_cast<size_t>(sample_y), color);
            }
        }
    }*/

    // extra parts
    // here I get the idea of step increment and design the frame. 
    // I write one row like float A01 = y0 - y1, B01 = x1 - x0; then AI help changed  the variable name and do repeated work
    /*float A01 = y0 - y1, B01 = x1 - x0;
    float A12 = y1 - y2, B12 = x2 - x1;
    float A20 = y2 - y0, B20 = x0 - x2;

    for(float x = leftx; x <= rightx; ++x) {
        for(float y = lowy; y <= highy; ++y) {
            float sample_x = x + 0.5f;
            float sample_y = y + 0.5f;

            float w0 = A12 * (sample_x - x2) + B12 * (sample_y - y2);
            float w1 = A20 * (sample_x - x0) + B20 * (sample_y - y0);
            float w2 = A01 * (sample_x - x1) + B01 * (sample_y - y1);

            if(w0 >= 0 && w1 >= 0 && w2 >= 0) {
                fill_pixel(static_cast<size_t>(x), static_cast<size_t>(y), color);
            }
        }
    }*/

    
            // TODO: Task 2: Update to implement super-sampled rasterization

    for (float x = leftx; x <= rightx; ++x) {
        for (float y = lowy; y <= highy; ++y) {
            
            // ed
            double rate_div = sqrt(sample_rate);
            /*float biass = 0.5f / rate_div;*/

            /*float sample_step = 1.0f / rate_div*/;

            float left_border_x = rate_div * x;
            float left_border_y = rate_div * y;

            float right_border_x = rate_div * (x + 1);
            float right_border_y = rate_div * (y + 1);

            for (float temp_x= left_border_x; temp_x < right_border_x; temp_x ++) {
                for (float temp_y = left_border_y; temp_y < right_border_y; temp_y ++) {
                    float sample_x = temp_x + 0.5f;
                    float sample_y = temp_y + 0.5f;

                    bool inside = ((x1 - x0) * (sample_y - rate_div * y0) - (y1 - y0) * (sample_x - rate_div * x0)) >= 0 && ((x2 - x1) * (sample_y - rate_div * y1) - (y2 - y1) * (sample_x - rate_div * x1)) >= 0 && ((x0 - x2) * (sample_y - rate_div * y2) - (y0 - y2) * (sample_x - rate_div *x2)) >= 0;

                    if (inside) {
                        fill_pixel(static_cast<size_t>(temp_x), static_cast<size_t>(temp_y), color);
                    }
                }
            }

            
        }
    }
  }


  void RasterizerImp::rasterize_interpolated_color_triangle(float x0, float y0, Color c0,
    float x1, float y1, Color c1,
    float x2, float y2, Color c2)
  {
    // TODO: Task 4: Rasterize the triangle, calculating barycentric coordinates and using them to interpolate vertex colors across the triangle
    // Hint: You can reuse code from rasterize_triangle
      float direction = (x1 - x0) * (y2 - y0) - (x2 - x0) * (y1 - y0);
      if (direction < 0) {
          swap(x1, x2);
          swap(y1, y2);
          swap(c1, c2);
      }

      // find many white line inside the image. refer to ed81a, use floor and ceil to solve
      float leftx = floor(min({x0, x1, x2}));
      float lowy = floor(min({y0, y1, y2}));

      float rightx = ceil(max({x0, x1, x2}));
      float highy = ceil(max({y0, y1, y2}));
      // ed
      double rate_div = sqrt(sample_rate);
      double area = abs((x1 - x0) * (y2 - y0) - (x2 - x0) * (y1 - y0));
      area *= rate_div;

      for (float x = leftx; x <= rightx; ++x) {
          for (float y = lowy; y <= highy; ++y) {

              
              /*float biass = 0.5f / rate_div;*/

              /*float sample_step = 1.0f / rate_div*/;

              float left_border_x = rate_div * x;
              float left_border_y = rate_div * y;

              float right_border_x = rate_div * (x + 1);
              float right_border_y = rate_div * (y + 1);

              for (float temp_x = left_border_x; temp_x < right_border_x; temp_x++) {
                  for (float temp_y = left_border_y; temp_y < right_border_y; temp_y++) {
                      float sample_x = temp_x + 0.5f;
                      float sample_y = temp_y + 0.5f;

                      bool inside = ((x1 - x0) * (sample_y - rate_div * y0) - (y1 - y0) * (sample_x - rate_div * x0)) >= 0 && ((x2 - x1) * (sample_y - rate_div * y1) - (y2 - y1) * (sample_x - rate_div * x1)) >= 0 && ((x0 - x2) * (sample_y - rate_div * y2) - (y0 - y2) * (sample_x - rate_div * x2)) >= 0;

                      if (inside) {

                          float w0 = ((y1 - y2) * (sample_x - rate_div *  x2) + (x2 - x1) * (sample_y - rate_div * y2)) / area;
                          float w1 = ((y2 - y0) * (sample_x - rate_div *  x2) + (x0 - x2) * (sample_y - rate_div * y2)) / area;
                          float w2 = 1 - w0 - w1;
                          Color actual_Color = w0 * c0 + w1 * c1 + w2 * c2;

                          fill_pixel(static_cast<size_t>(temp_x), static_cast<size_t>(temp_y), actual_Color);
                      }
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
      /*int level = tex.get_level;*/
      int level = 0;
      float direction = (x1 - x0) * (y2 - y0) - (x2 - x0) * (y1 - y0);
      if (direction < 0) {
          swap(x1, x2);
          swap(y1, y2);
          swap(u1, u2); 
          swap(v1, v2);
      }

      // find many white line inside the image. refer to ed81a, use floor and ceil to solve
      float leftx = floor(min({x0, x1, x2}));
      float lowy = floor(min({y0, y1, y2}));

      float rightx = ceil(max({x0, x1, x2}));
      float highy = ceil(max({y0, y1, y2}));
      // ed
      double rate_div = sqrt(sample_rate);
      double area = abs((x1 - x0) * (y2 - y0) - (x2 - x0) * (y1 - y0));
      area *= rate_div;

      for(float x = leftx; x <= rightx; ++x) {
          for(float y = lowy; y <= highy; ++y) {


              float left_border_x = rate_div * x;
              float left_border_y = rate_div * y;

              float right_border_x = rate_div * (x + 1);
              float right_border_y = rate_div * (y + 1);

              for(float temp_x = left_border_x; temp_x < right_border_x; temp_x++) {
                  for(float temp_y = left_border_y; temp_y < right_border_y; temp_y++) {
                      float sample_x = temp_x + 0.5f;
                      float sample_y = temp_y + 0.5f;

                      bool inside = ((x1 - x0) * (sample_y - rate_div * y0) - (y1 - y0) * (sample_x - rate_div * x0)) >= 0 && ((x2 - x1) * (sample_y - rate_div * y1) - (y2 - y1) * (sample_x - rate_div * x1)) >= 0 && ((x0 - x2) * (sample_y - rate_div * y2) - (y0 - y2) * (sample_x - rate_div * x2)) >= 0;

                      if(inside) {

                          //calculate the x,y
                          float w0 = ((y1 - y2) * (sample_x - rate_div * x2) + (x2 - x1) * (sample_y - rate_div * y2)) / area;
                          float w1 = ((y2 - y0) * (sample_x - rate_div * x2) + (x0 - x2) * (sample_y - rate_div * y2)) / area;
                          float w2 = 1 - w0 - w1;

                          // u and v
                          float u = w0 * u0 + w1 * u1 + w2 * u2;
                          float v = w0 * v0 + w1 * v1 + w2 * v2;

                          // x+1,y
                          float sample_x_plus_1 = temp_x + 1.0f + 0.5f; 
                          float sample_y_same = temp_y + 0.5f; 

                          float w0_x_plus_1 = ((y1 - y2) * (sample_x_plus_1 - rate_div * x2) + (x2 - x1) * (sample_y_same - rate_div * y2)) / area;
                          float w1_x_plus_1 = ((y2 - y0) * (sample_x_plus_1 - rate_div * x2) + (x0 - x2) * (sample_y_same - rate_div * y2)) / area;
                          float w2_x_plus_1 = 1 - w0_x_plus_1 - w1_x_plus_1;

                          // u'
                          float u_x_plus_1 = w0_x_plus_1 * u0 + w1_x_plus_1 * u1 + w2_x_plus_1 * u2;

                          bool inside_x_plus_1 = w0_x_plus_1 >= 0 && w1_x_plus_1 >= 0 && w2_x_plus_1 >= 0;
                          if (!inside_x_plus_1) { 
                              sample_x_plus_1 = temp_x - 1.0f + 0.5f; 
                               // x-1
                              w0_x_plus_1 = ((y1 - y2) * (sample_x_plus_1 - rate_div * x2) + (x2 - x1) * (sample_y_same - rate_div * y2)) / area;
                              w1_x_plus_1 = ((y2 - y0) * (sample_x_plus_1 - rate_div * x2) + (x0 - x2) * (sample_y_same - rate_div * y2)) / area;
                              w2_x_plus_1 = 1 - w0_x_plus_1 - w1_x_plus_1;
                              u_x_plus_1 = w0_x_plus_1 * u0 + w1_x_plus_1 * u1 + w2_x_plus_1 * u2; // ¸üÐÂu_x_plus_1
                          }

                          // x,y+1
                          float sample_x_same = temp_x + 0.5f;
                          float sample_y_plus_1 = temp_y + 1.0f + 0.5f; 


                          float w0_y_plus_1 = ((y1 - y2) * (sample_x_same - rate_div * x2) + (x2 - x1) * (sample_y_plus_1 - rate_div * y2)) / area;
                          float w1_y_plus_1 = ((y2 - y0) * (sample_x_same - rate_div * x2) + (x0 - x2) * (sample_y_plus_1 - rate_div * y2)) / area;
                          float w2_y_plus_1 = 1 - w0_y_plus_1 - w1_y_plus_1;

                          // v'
                          float v_y_plus_1 = w0_y_plus_1 * v0 + w1_y_plus_1 * v1 + w2_y_plus_1 * v2;

                          bool inside_y_plus_1 = w0_y_plus_1 >= 0 && w1_y_plus_1 >= 0 && w2_y_plus_1 >= 0;
                          if(!inside_y_plus_1) { 
                              sample_y_plus_1 = temp_y - 1.0f + 0.5f; 
                              // y-1
                              w0_y_plus_1 = ((y1 - y2) * (sample_x_same - rate_div * x2) + (x2 - x1) * (sample_y_plus_1 - rate_div * y2)) / area;
                              w1_y_plus_1 = ((y2 - y0) * (sample_x_same - rate_div * x2) + (x0 - x2) * (sample_y_plus_1 - rate_div * y2)) / area;
                              w2_y_plus_1 = 1 - w0_y_plus_1 - w1_y_plus_1;
                              v_y_plus_1 = w0_y_plus_1 * v0 + w1_y_plus_1 * v1 + w2_y_plus_1 * v2;
                          }

                          float ddx_u = u_x_plus_1 - u; 
                          float ddy_v = v_y_plus_1 - v; 

                          SampleParams sp;
                          sp.lsm = lsm;
                          sp.psm = psm;

                          sp.p_uv = Vector2D(u, v);

                          sp.p_dx_uv = Vector2D(ddx_u, 0); 
                          sp.p_dy_uv = Vector2D(0, ddy_v); 
                          

                          
                          Color actual_Color = tex.sample(sp);

                          

                          fill_pixel(static_cast<size_t>(temp_x), static_cast<size_t>(temp_y), actual_Color);
                      }
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
    this->sample_buffer.resize(width * height * rate, Color::White);
    // clear buffer
    clear_buffers();
  }


  void RasterizerImp::set_framebuffer_target(unsigned char* rgb_framebuffer,
    size_t width, size_t height)
  {
    // TODO: Task 2: You may want to update this function for supersampling support
    

    this->width = width;
    this->height = height;
    this->rgb_framebuffer_target = rgb_framebuffer;

    

    this->sample_buffer.resize(width * height * sample_rate, Color::White);

    // clear buffer
    clear_buffers();
  }


  void RasterizerImp::clear_buffers() {
    std::fill(rgb_framebuffer_target, rgb_framebuffer_target + 3 * width * height, 255);
    std::fill(sample_buffer.begin(), sample_buffer.end(), Color::White);
  }


  // This function is called at the end of rasterizing all elements of the
  // SVG file.  If you use a supersample buffer to rasterize SVG elements
  // for antialising, you could use this call to fill the target framebuffer
  // pixels from the supersample buffer data.
  //
  void RasterizerImp::resolve_to_framebuffer() {
    
      // TODO: Task 2: You will likely want to update this function for supersampling support



    /*for (int x = 0; x < width; ++x) {
      for (int y = 0; y < height; ++y) {
        Color col = sample_buffer[y * width + x];

        for (int k = 0; k < 3; ++k) {
          this->rgb_framebuffer_target[3 * (y * width + x) + k] = (&col.r)[k] * 255;
        }
      }
    }*/
      double rate_div = sqrt(sample_rate);
      for (size_t x = 0; x < width; ++x) {
          for (size_t y = 0; y < height; ++y) {

              // ed
              
              /*float biass = 0.5f / rate_div;*/

              /*float sample_step = 1.0f / rate_div*/;

              float left_border_x = rate_div * x;
              float left_border_y = rate_div * y;

              float right_border_x = rate_div * (x + 1);
              float right_border_y = rate_div * (y + 1);

              Color avg_color(0, 0, 0);
              for (float temp_x = left_border_x; temp_x < right_border_x; temp_x++) {
                  for (float temp_y = left_border_y; temp_y < right_border_y; temp_y++) {
                      size_t sample_index = (temp_y * width * rate_div) + temp_x;
                      avg_color += sample_buffer[sample_index];
                  }
              }

              avg_color.r /= sample_rate;
              avg_color.g /= sample_rate;
              avg_color.b /= sample_rate;
              for (int k = 0; k < 3; ++k) {
                  this->rgb_framebuffer_target[3 * (y * width + x) + k] = (&avg_color.r)[k] * 255;
              }

          }
      }

  }

  Rasterizer::~Rasterizer() { }


}// CGL
