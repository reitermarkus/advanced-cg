mod vector;
mod triangular;
mod triangle;
mod image;
mod patch_triangle;
mod ray;

use std::vec::Vec;
use std::collections::HashMap;
use std::f64::consts::PI;
use std::io::{self, Write};

use image::Image;
use ray::Ray;
use triangular::Triangular;
use triangle::Triangle;
use patch_triangle::PatchTriangle;
use vector::Vector;

fn tris() -> Vec<Triangle> {
  vec![
  /* Cornell Box walls */
    Triangle::new(Vector::new(  0.0,  0.0,   0.0), Vector::new( 100.0, 0.0,    0.0), Vector::new(0.0,  80.0,    0.0), Vector::new(0.0, 0.0, 0.0), Vector::new(0.75, 0.75, 0.75)), // Back:   bottom-left
    Triangle::new(Vector::new(100.0, 80.0,   0.0), Vector::new(-100.0, 0.0,    0.0), Vector::new(0.0, -80.0,    0.0), Vector::new(0.0, 0.0, 0.0), Vector::new(0.75, 0.75, 0.75)), // Back:   top-right
    Triangle::new(Vector::new(  0.0,  0.0, 170.0), Vector::new( 100.0, 0.0,    0.0), Vector::new(0.0,   0.0, -170.0), Vector::new(0.0, 0.0, 0.0), Vector::new(0.75, 0.75, 0.75)), // Bottom: front-left
    Triangle::new(Vector::new(100.0,  0.0,   0.0), Vector::new(-100.0, 0.0,    0.0), Vector::new(0.0,   0.0,  170.0), Vector::new(0.0, 0.0, 0.0), Vector::new(0.75, 0.75, 0.75)), // Bottom: back-right
    Triangle::new(Vector::new(  0.0, 80.0,   0.0), Vector::new( 100.0, 0.0,    0.0), Vector::new(0.0,   0.0,  170.0), Vector::new(0.0, 0.0, 0.0), Vector::new(0.75, 0.75, 0.75)), // Top:    back-left
    Triangle::new(Vector::new(100.0, 80.0, 170.0), Vector::new(-100.0, 0.0,    0.0), Vector::new(0.0,   0.0, -170.0), Vector::new(0.0, 0.0, 0.0), Vector::new(0.75, 0.75, 0.75)), // Top:    front-right
    Triangle::new(Vector::new(  0.0,  0.0, 170.0), Vector::new(   0.0, 0.0, -170.0), Vector::new(0.0,  80.0,    0.0), Vector::new(0.0, 0.0, 0.0), Vector::new(0.75, 0.25, 0.25)), // Left:   front-bottom
    Triangle::new(Vector::new(  0.0, 80.0,   0.0), Vector::new(   0.0, 0.0,  170.0), Vector::new(0.0, -80.0,    0.0), Vector::new(0.0, 0.0, 0.0), Vector::new(0.75, 0.25, 0.25)), // Left:   back-top
    Triangle::new(Vector::new(100.0,  0.0,   0.0), Vector::new(   0.0, 0.0,  170.0), Vector::new(0.0,  80.0,    0.0), Vector::new(0.0, 0.0, 0.0), Vector::new(0.25, 0.25, 0.75)), // Right:  back-bottom
    Triangle::new(Vector::new(100.0, 80.0, 170.0), Vector::new(   0.0, 0.0, -170.0), Vector::new(0.0, -80.0,    0.0), Vector::new(0.0, 0.0, 0.0), Vector::new(0.25, 0.25, 0.75)), // Right:  front-top
    Triangle::new(Vector::new(100.0,  0.0, 170.0), Vector::new(-100.0, 0.0,    0.0), Vector::new(0.0,  80.0,    0.0), Vector::new(0.0, 0.0, 0.0), Vector::new(0.0,  1.0,  0.0)),  // Front:  bottom-right (not visible)
    Triangle::new(Vector::new(  0.0, 80.0, 170.0), Vector::new( 100.0, 0.0,    0.0), Vector::new(0.0, -80.0,    0.0), Vector::new(0.0, 0.0, 0.0), Vector::new(0.0,  1.0,  0.0)),  // Front:  top-left (not visible)

    /* Area light source on top */
    Triangle::new(Vector::new(40.0, 79.99, 65.0), Vector::new( 20.0, 0.0, 0.0), Vector::new(0.0, 0.0,  20.0), Vector::new(12.0, 12.0, 12.0), Vector::new(0.75, 0.75, 0.75)), // back-left
    Triangle::new(Vector::new(60.0, 79.99, 85.0), Vector::new(-20.0, 0.0, 0.0), Vector::new(0.0, 0.0, -20.0), Vector::new(12.0, 12.0, 12.0), Vector::new(0.75, 0.75, 0.75)), // front-right

    /* Cuboid in room */
    Triangle::new(Vector::new(30.0,  0.0, 100.0), Vector::new(  0.0, 0.0, -20.0), Vector::new(0.0,  40.0,   0.0), Vector::new(0.0, 0.0, 0.0), Vector::new(0.75, 0.75, 0.75)), // Right: front-bottom
    Triangle::new(Vector::new(30.0, 40.0,  80.0), Vector::new(  0.0, 0.0,  20.0), Vector::new(0.0, -40.0,   0.0), Vector::new(0.0, 0.0, 0.0), Vector::new(0.75, 0.75, 0.75)), // Right: back-top
    Triangle::new(Vector::new(10.0,  0.0,  80.0), Vector::new(  0.0, 0.0,  20.0), Vector::new(0.0,  40.0,   0.0), Vector::new(0.0, 0.0, 0.0), Vector::new(0.75, 0.75, 0.75)), // Left:  back-bottom
    Triangle::new(Vector::new(10.0, 40.0, 100.0), Vector::new(  0.0, 0.0, -20.0), Vector::new(0.0, -40.0,   0.0), Vector::new(0.0, 0.0, 0.0), Vector::new(0.75, 0.75, 0.75)), // Left:  front-top
    Triangle::new(Vector::new(10.0,  0.0, 100.0), Vector::new( 20.0, 0.0,   0.0), Vector::new(0.0,  40.0,   0.0), Vector::new(0.0, 0.0, 0.0), Vector::new(0.75, 0.75, 0.75)), // Front: bottom-left
    Triangle::new(Vector::new(30.0, 40.0, 100.0), Vector::new(-20.0, 0.0,   0.0), Vector::new(0.0, -40.0,   0.0), Vector::new(0.0, 0.0, 0.0), Vector::new(0.75, 0.75, 0.75)), // Front: top-right
    Triangle::new(Vector::new(30.0,  0.0,  80.0), Vector::new(-20.0, 0.0,   0.0), Vector::new(0.0,  40.0,   0.0), Vector::new(0.0, 0.0, 0.0), Vector::new(0.75, 0.75, 0.75)), // Back:  bottom-right
    Triangle::new(Vector::new(10.0,  4.0,  80.0), Vector::new( 20.0, 0.0,   0.0), Vector::new(0.0, -40.0,   0.0), Vector::new(0.0, 0.0, 0.0), Vector::new(0.75, 0.75, 0.75)), // Back:  top-left
    Triangle::new(Vector::new(10.0, 40.0, 100.0), Vector::new( 20.0, 0.0,   0.0), Vector::new(0.0,   0.0, -20.0), Vector::new(0.0, 0.0, 0.0), Vector::new(0.75, 0.75, 0.75)), // Top:   front-left
    Triangle::new(Vector::new(30.0, 40.0,  80.0), Vector::new(-20.0, 0.0,   0.0), Vector::new(0.0,   0.0,  20.0), Vector::new(0.0, 0.0, 0.0), Vector::new(0.75, 0.75, 0.75)), // Top:   back-right
  ]
}

fn intersect_scene(tris: &[Triangle], ray: &Ray, t: &mut f64, id: &mut i64, normal: &mut Vector) -> bool {
  *t = 1e20;
  *id = -1;

  let n = tris.len();
  for i in 0..n {
    let d = tris[i].intersect(ray);
    if d > 0.0 && d < *t {
      *t = d;
      *id = i as i64;
      *normal = tris[i].normal;
    }
  }

  *t < 1e20
}

fn calculate_form_factors(tris: &mut [Triangle], divisions: u64, mc_sample: i64) -> HashMap<usize, Vec<HashMap<usize, Vec<f64>>>> {
  let mut form_factors: HashMap<usize, Vec<HashMap<usize, Vec<f64>>>> = HashMap::new();
  let mut patch_num = 0;

  // Total number of patches in scene.
  let n = tris.len();
  for i in 0..n {
    tris[i].init_patches(divisions);
    patch_num += tris[i].patches.len();
  }

  let form_factor_num = patch_num.pow(2);

  println!("Number of triangles: {}", n);
  println!("Number of patches: {}", patch_num);
  println!("Number of form factors: {}", form_factor_num);

  for i in 0..n {
    form_factors.insert(i, vec![HashMap::new(); tris[i].patches.len()]);

    for p in 0..tris[i].patches.len() {
      for j in 0..n {
        form_factors.get_mut(&i).unwrap()[p].insert(j, vec![0.0; tris[j].patches.len()]);
      }
    }
  }

  let mut patch_area: HashMap<usize, Vec<f64>> = HashMap::new();

  // Precompute patch areas, assuming same size for each triangle.
  for i in 0..n {
    let area = tris[i].area / tris[i].patches.len() as f64;
    patch_area.insert(i, vec![area; tris[i].patches.len()]);
  }

  // Loop over all triangles in scene.
  for i in 0..n {
    print!("{} ", i);

    // Loop over all patches in rectangle i.
    for p_i in 0..tris[i].patches.len() {
      if p_i % tris[i].divisions as usize == 0 {
        print!("*");
        io::stdout().flush().unwrap();
      }

      // Loop over all triangles in scene for triangles i.
      for j in 0..n {
        // Loop over all patches in rectangle j.
        for p_j in 0..tris[j].patches.len() {
          // Do not compute form factors for patches on same rectangle;
          // also exploit symmetry to reduce computation;
          // intemediate values; will be divided by patch area below.
          if i < j {
            let mut form_factor = 0.0;

            // Monte Carlo integration of form factor double integral.

            // Uniform PDF for Monte Carlo (1 / Ai) x (1 / Aj).
            let pdf = (1.0 / patch_area[&i][p_i]) *
                      (1.0 / patch_area[&j][p_j]);

            // Determine rays of NixNi uniform samples of patch
            // on i to NjxNj uniform samples of patch on j.
            for ii in 0..mc_sample {
              for jj in 0..mc_sample {
                let t_i: &PatchTriangle = &tris[i].sub_triangles[p_i];
                let t_j: &PatchTriangle = &tris[j].sub_triangles[p_j];

                let xi: Vector = t_i.random_sample();
                let xj: Vector = t_j.random_sample();

                // Check for visibility between sample points.
                let ij: Vector = (&xj - &xi).normalize();

                let mut t = 0.0;
                let mut id = -1;
                let mut normal = Vector::new(0.0, 0.0, 0.0);
                if intersect_scene(&tris, &Ray::new(&xi, &ij), &mut t, &mut id, &mut normal) && id != j as i64 {
                  continue; // If intersection with other triangle.
                }

                // Cosines of angles beteen normals and ray inbetween.
                let d0 = tris[i].normal.dot_product(&ij);
                let d1 = tris[j].normal.dot_product(&(-1.0 * ij));

                // Continue if patches facing each other.
                if d0 > 0.0 && d1 > 0.0 {
                  // Sample form factor.
                  let k = d0 * d1 / (PI * (&xj - &xi).length_squared());

                  // Add weighted sample to estimate.
                  form_factor += k / pdf;
                }
              }
            }

            // Divide by number of samples.
            form_factor /= mc_sample.pow(2) as f64;

            form_factors.get_mut(&i).unwrap()[p_i].get_mut(&j).unwrap().insert(p_j, form_factor);
            form_factors.get_mut(&j).unwrap()[p_j].get_mut(&i).unwrap().insert(p_i, form_factor);
          }
        }
      }
    }

    println!("");
  }

  // Divide by area to get final form factors.
  for i in 0..n {
    for p_a in 0..tris[i].patches.len() {
      for j in 0..n {
        for p_b in 0..tris[j].patches.len() {
          let area = patch_area[&i][p_a];
          let ff = form_factors[&i][p_a][&j][p_b] / area;
          form_factors.get_mut(&i).unwrap()[p_a].get_mut(&j).unwrap().insert(p_b, if ff < 0.0 { 0.0 } else if ff > 1.0 { 1.0 } else { ff });
        }
      }
    }
  }

  form_factors
}

fn calculate_radiosity(tris: &mut [Triangle], form_factors: &HashMap<usize, Vec<HashMap<usize, Vec<f64>>>>) {
  for i in 0..tris.len() {
    for p_a in 0..tris[i].patches.len() {
      let mut color = Vector::new(0.0, 0.0, 0.0);

      for j in 0..tris.len() {
        for p_b in 0..tris[j].patches.len() {
          let form_factor_ij = form_factors[&i][p_a][&j][p_b];

          // Add form factor multiplied with radiosity of previous step.
          if form_factor_ij > 0.0 {
            color = color + form_factor_ij * tris[j].patches[p_b];
          }
        }
      }

      // Multiply sum with color of patch and add emission.
      color = tris[i].color.entrywise_product(&color) + tris[i].emission;

      // Store overall patch radiosity of current iteration.
      tris[i].patches[p_a] = color;
    }
  }
}

fn main() {
  let mut tris = tris();

  let divisions: u64 = 8;
  let samples = 4;

  let form_factors = calculate_form_factors(&mut tris, divisions, samples);

  calculate_radiosity(&mut tris, &form_factors);

  let image_width = 640;
  let image_height = 480;

  let image = Image::new(image_width, image_height);


  let v1 = Vector::new(1.0, 2.0, 3.0);
  let v2 = Vector::new(3.0, 2.0, 1.0);


  println!("Hello, world!");
  println!("{:?}", (v1 + v2));


  if let Err(e) = image.save(&"yolo.ppm".to_string()) {
    panic!("{:?}", e)
  }
}
