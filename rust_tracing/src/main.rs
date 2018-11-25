use std::f64::consts::PI;
use std::thread;
use std::sync::mpsc::channel;

#[macro_use]
extern crate clap;
use clap::{App, Arg};
extern crate crossbeam;
#[macro_use]
extern crate downcast_rs;
#[macro_use]
extern crate lazy_static;
extern crate rand;
extern crate rayon;
use rayon::prelude::*;

extern crate image;

extern crate indicatif;
use indicatif::{ProgressBar, ProgressStyle};

#[macro_use]
extern crate glium;

use glium::{
  glutin::{Api, GlProfile, GlRequest}, index::{NoIndices, PrimitiveType},
  texture::buffer_texture::{BufferTexture, BufferTextureType}, vertex::EmptyVertexAttributes,
  Surface,
};

#[macro_use]
extern crate glium_text_rusttype as glium_text;

mod color;
use color::Color;
mod ray;
use ray::Ray;
mod sampler;
use sampler::{drand48, random_direction};
mod scene_object;
use scene_object::{ReflType, SceneObject};
mod sphere;
use sphere::Sphere;
mod triangle;
use triangle::Triangle;
mod vector;
use vector::Vector;

lazy_static! {
  static ref TRIS: Vec<Triangle> = {
    vec![
      /* Cornell Box walls */
      Triangle::new(Vector::new(  0.0,  0.0,   0.0), Vector::new( 100.0, 0.0,    0.0), Vector::new(0.0,  80.0,    0.0), Color::zero(), Color::new(0.75, 0.75, 0.75), ReflType::DIFF), // Back:   bottom-left
      Triangle::new(Vector::new(100.0, 80.0,   0.0), Vector::new(-100.0, 0.0,    0.0), Vector::new(0.0, -80.0,    0.0), Color::zero(), Color::new(0.75, 0.75, 0.75), ReflType::DIFF), // Back:   top-right
      Triangle::new(Vector::new(  0.0,  0.0, 170.0), Vector::new( 100.0, 0.0,    0.0), Vector::new(0.0,   0.0, -170.0), Color::zero(), Color::new(0.75, 0.75, 0.75), ReflType::DIFF), // Bottom: front-left
      Triangle::new(Vector::new(100.0,  0.0,   0.0), Vector::new(-100.0, 0.0,    0.0), Vector::new(0.0,   0.0,  170.0), Color::zero(), Color::new(0.75, 0.75, 0.75), ReflType::DIFF), // Bottom: back-right
      Triangle::new(Vector::new(  0.0, 80.0,   0.0), Vector::new( 100.0, 0.0,    0.0), Vector::new(0.0,   0.0,  170.0), Color::zero(), Color::new(0.75, 0.75, 0.75), ReflType::DIFF), // Top:    back-left
      Triangle::new(Vector::new(100.0, 80.0, 170.0), Vector::new(-100.0, 0.0,    0.0), Vector::new(0.0,   0.0, -170.0), Color::zero(), Color::new(0.75, 0.75, 0.75), ReflType::DIFF), // Top:    front-right
      Triangle::new(Vector::new(  0.0,  0.0, 170.0), Vector::new(   0.0, 0.0, -170.0), Vector::new(0.0,  80.0,    0.0), Color::zero(), Color::new(0.75, 0.25, 0.25), ReflType::DIFF), // Left:   front-bottom
      Triangle::new(Vector::new(  0.0, 80.0,   0.0), Vector::new(   0.0, 0.0,  170.0), Vector::new(0.0, -80.0,    0.0), Color::zero(), Color::new(0.75, 0.25, 0.25), ReflType::DIFF), // Left:   back-top
      Triangle::new(Vector::new(100.0,  0.0,   0.0), Vector::new(   0.0, 0.0,  170.0), Vector::new(0.0,  80.0,    0.0), Color::zero(), Color::new(0.25, 0.25, 0.75), ReflType::DIFF), // Right:  back-bottom
      Triangle::new(Vector::new(100.0, 80.0, 170.0), Vector::new(   0.0, 0.0, -170.0), Vector::new(0.0, -80.0,    0.0), Color::zero(), Color::new(0.25, 0.25, 0.75), ReflType::DIFF), // Right:  front-top
      Triangle::new(Vector::new(100.0,  0.0, 170.0), Vector::new(-100.0, 0.0,    0.0), Vector::new(0.0,  80.0,    0.0), Color::zero(), Color::new(0.25, 0.75,  0.25), ReflType::DIFF),  // Front:  bottom-right (not visible)
      Triangle::new(Vector::new(  0.0, 80.0, 170.0), Vector::new( 100.0, 0.0,    0.0), Vector::new(0.0, -80.0,    0.0), Color::zero(), Color::new(0.25, 0.75,  0.25), ReflType::DIFF),  // Front:  top-left (not visible)

      /* Cuboid in room */
      Triangle::new(Vector::new(30.0,  0.0, 100.0), Vector::new(  0.0, 0.0, -20.0), Vector::new(0.0,  40.0,   0.0), Color::zero(), Color::new(0.75, 0.75, 0.75), ReflType::TRAN), // Right: front-bottom
      Triangle::new(Vector::new(30.0, 40.0,  80.0), Vector::new(  0.0, 0.0,  20.0), Vector::new(0.0, -40.0,   0.0), Color::zero(), Color::new(0.75, 0.75, 0.75), ReflType::TRAN), // Right: back-top
      Triangle::new(Vector::new(10.0,  0.0,  80.0), Vector::new(  0.0, 0.0,  20.0), Vector::new(0.0,  40.0,   0.0), Color::zero(), Color::new(0.75, 0.75, 0.75), ReflType::TRAN), // Left:  back-bottom
      Triangle::new(Vector::new(10.0, 40.0, 100.0), Vector::new(  0.0, 0.0, -20.0), Vector::new(0.0, -40.0,   0.0), Color::zero(), Color::new(0.75, 0.75, 0.75), ReflType::TRAN), // Left:  front-top
      Triangle::new(Vector::new(10.0,  0.0, 100.0), Vector::new( 20.0, 0.0,   0.0), Vector::new(0.0,  40.0,   0.0), Color::zero(), Color::new(0.75, 0.75, 0.75), ReflType::TRAN), // Front: bottom-left
      Triangle::new(Vector::new(30.0, 40.0, 100.0), Vector::new(-20.0, 0.0,   0.0), Vector::new(0.0, -40.0,   0.0), Color::zero(), Color::new(0.75, 0.75, 0.75), ReflType::TRAN), // Front: top-right
      Triangle::new(Vector::new(30.0,  0.0,  80.0), Vector::new(-20.0, 0.0,   0.0), Vector::new(0.0,  40.0,   0.0), Color::zero(), Color::new(0.75, 0.75, 0.75), ReflType::TRAN), // Back:  bottom-right
      Triangle::new(Vector::new(10.0, 40.0,  80.0), Vector::new( 20.0, 0.0,   0.0), Vector::new(0.0, -40.0,   0.0), Color::zero(), Color::new(0.75, 0.75, 0.75), ReflType::TRAN), // Back:  top-left
      Triangle::new(Vector::new(10.0, 40.0, 100.0), Vector::new( 20.0, 0.0,   0.0), Vector::new(0.0,   0.0, -20.0), Color::zero(), Color::new(0.75, 0.75, 0.75), ReflType::TRAN), // Top:   front-left
      Triangle::new(Vector::new(30.0, 40.0,  80.0), Vector::new(-20.0, 0.0,   0.0), Vector::new(0.0,   0.0,  20.0), Color::zero(), Color::new(0.75, 0.75, 0.75), ReflType::TRAN), // Top:   back-right
    ]
  };

  static ref SPHERES: Vec<Sphere> = {
    vec![
      Sphere::new(16.5, Vector::new(27.0, 16.5, 47.0), Color::zero(), Color::new(1.0, 1.0, 1.0), ReflType::SPEC),  // Mirror Sphere
      Sphere::new(16.5, Vector::new(73.0, 16.5, 78.0), Color::zero(), Color::new(1.0, 1.0, 1.0), ReflType::REFR),  // Mirror Sphere
      Sphere::new(12.5, Vector::new(48.0, 12.5, 117.0), Color::zero(), Color::new(1.0, 1.0, 1.0), ReflType::GLOS),  // Mirror Sphere

      Sphere::new(1.5, Vector::new(50.0, 81.6 - 16.5, 81.6), Color::new(4.0, 4.0, 4.0) * 100.0, Color::zero(), ReflType::DIFF) // Light
    ]
  };
}

fn intersect_scene(objects: &[Box<&dyn SceneObject>], ray: &Ray, t: &mut f64, id: &mut usize) -> bool {
  *t = 1e20;

  for i in 0..objects.len() {
    let d = objects[i].intersect(ray);
    if d > 0.0 && d < *t {
      *t = d;
      *id = i;
    }
  }

  *t < 1e20
}

fn perfect_reflection(dir: &Vector, normal: &Vector) -> Vector {
  *dir - *normal * 2.0 * (*normal).dot_product(dir)
}

fn radiance(objects: &[Box<&dyn SceneObject>], ray: &Ray, mut depth: i32, q: i32) -> Color {
  depth += 1;
  let mut t = 0.0;
  let mut id = 0;

  if !intersect_scene(objects, ray, &mut t, &mut id) {
     /* No intersection with scene */
    return Color::new(0.0, 0.0, 0.0);
  }

  /* Intersection point */
  let hitpoint = ray.org + ray.dir * t;
  let obj = &objects[id];

  /* Normal at intersection */
  let normal;

  if let Some(sphere) = obj.downcast_ref::<Sphere>() {
    normal = (hitpoint - sphere.position).normalize();
  } else if let Some(triangle) = obj.downcast_ref::<Triangle>() {
    normal = triangle.normal;
  } else {
    panic!("Object is not a sphere nor a triangle!")
  }

  let mut nl = normal;

  /* Obtain flipped normal, if object hit from inside */
  if normal.dot_product(&ray.dir) >= 0.0 {
    nl = -nl;
  }

  let mut col = obj.color();

  /* Maximum RGB reflectivity for Russian Roulette */
  let p = col.clamp(0.0, 0.999).max();

  /* After 5 bounces or if max reflectivity is zero */
  if depth > 5 || p == 0.0 {
    if drand48() >= p {
      /* No further bounces, only return potential emission */
      return obj.emission() * q as f64;
    }

    /* Scale estimator to remain unbiased */
    col = col * (1.0 / p);
  }

  let mut ray_dir = ray.dir;
  let mut nt = 1.5;

  match obj.refl() {
    ReflType::DIFF => {
      /* Compute random reflection vector on hemisphere */
      let r1 = 2.0 * PI * drand48();
      let r2 = drand48();
      let r2s = r2.sqrt();

      /* Set up local orthogonal coordinate system u,v,w on surface */
      let w = nl;
      let u = if w.x.abs() > 0.1 { Vector::new(0.0, 1.0, 0.0) } else { Vector::new(1.0, 0.0, 0.0) };
      let u = u.cross_product(&w).normalize();

      let v = w.cross_product(&u);

      /* Random reflection vector d */
      let d = (u * r1.cos() * r2s +
               v * r1.sin() * r2s +
               w * (1.0 - r2).sqrt()).normalize();

      let mut e: Vector = Vector::zero();
      for i in 0..objects.len() {
        let light_source = &objects[i];

        if light_source.emission().x <= 0.0 && light_source.emission().y <= 0.0 && light_source.emission().z <= 0.0 {
          continue;
        }

        let sphere: &Sphere = match light_source.downcast_ref::<Sphere>() {
          Some(s) => s,
          None => {
            println!("Warning: Only spherical light sources are implemented.");
            continue;
          }
        };

        /* Randomly sample spherical light source from surface intersection */

        // Create random sample direction towards spherical light source.
        let cos_a_max = (1.0 - sphere.radius.powf(2.0) /
                                (hitpoint - sphere.position).dot_product(&(hitpoint - sphere.position))).sqrt();

        let l = random_direction(sphere.position - hitpoint, cos_a_max.acos());

        if intersect_scene(objects, &Ray::new(hitpoint, l), &mut t, &mut id) && id == i {
          let omega = 2.0 * PI * (1.0 - cos_a_max);
          e += col.entrywise_product(sphere.emission() * l.dot_product(&nl) * omega) / PI;
        }
      }

          /* Return potential light emission, direct lighting, and indirect lighting (via
        recursive call for Monte-Carlo integration */
      return obj.emission() * q as f64 + e + col.entrywise_product(radiance(objects, &Ray::new(hitpoint, d), depth, 0));
    },
    ReflType::GLOS => {
      let angle = 0.10;
      let l = random_direction(nl, angle);

      return obj.emission() + col.entrywise_product(radiance(objects, &Ray::new(hitpoint, l), depth, 1));
    },
    ReflType::SPEC => {
      // Return light emission mirror reflection (via recursive call using perfect reflection vector).
      let reflection_ray = Ray::new(hitpoint, perfect_reflection(&ray.dir, &normal));
      return obj.emission() + col.entrywise_product(radiance(objects, &reflection_ray, depth, 1));
    },
    ReflType::TRAN => {
      let angle = 0.2;
      ray_dir = random_direction(ray.dir, angle);
      nt = 1.15;
    },
    _ => {}
  }

  let nc = 1.0;

  let into = normal.dot_product(&nl) > 0.0;       /* Bool for checking if ray from outside going in */
  let nnt = if into { nc / nt } else { nt / nc }; /* Set ratio depending on hit from inside or outside */

  let ddn = ray_dir.dot_product(&nl);
  let cos2t = 1.0 - nnt * nnt * (1.0 - ddn * ddn);

  let reflection_ray = Ray::new(hitpoint, perfect_reflection(&ray_dir, &normal)); // Perfect reflection.

    // Check for total internal reflection, if so only reflect.
  if cos2t < 0.0 {
    return obj.emission() +  col.entrywise_product(radiance(objects, &reflection_ray, depth, 1));
  }

  // Determine transmitted ray direction for refraction.
  let transmission_direction = if into {
    ray_dir * nnt - normal * (ddn * nnt + cos2t.sqrt())
  } else {
    ray_dir * nnt + normal * (ddn * nnt + cos2t.sqrt())
  };

  let transmission_direction = transmission_direction.normalize();

  /* Cosine of correct angle depending on outside/inside */
  let c = if into { 1.0 + ddn } else { 1.0 - transmission_direction.dot_product(&normal) };

  /* Compute Schlick's approximation of Fresnel equation */
  let r0 = ((nc - nt) / (nc + nt)).powf(2.0);
  let re = r0 + (1.0 - r0) * c.powf(5.0);   // Reflectance
  let tr = 1.0 - re;                     // Transmittance

  let transmission_ray = Ray::new(hitpoint, transmission_direction);

  // Initially, use both reflection and trasmission.
  if depth < 3 {
    return obj.emission() +  col.entrywise_product(radiance(objects, &reflection_ray, depth, 1) * re
                                                   + radiance(objects, &transmission_ray, depth, 1) * tr);
  }
  // Probability for selecting reflectance or transmittance */
  let p = 0.25 + 0.5 * re;
  let rp = re / p;         /* Scaling factors for unbiased estimator */
  let tp = tr / (1.0 - p);

  /* Russian Roulette */
  if drand48() < p {
    return obj.emission() +  col.entrywise_product(radiance(objects, &reflection_ray, depth, 1) * rp);
  }

  return obj.emission() +  col.entrywise_product(radiance(objects, &transmission_ray, depth, 1) * tp);
}

fn main() {
  let matches = App::new("rust_tracing")
                  .arg(Arg::with_name("height")
                    .long("height")
                    .takes_value(true)
                    .min_values(1)
                    .max_values(1)
                    .multiple(true)
                    .help("Sets height of image"))
                  .arg(Arg::with_name("width")
                    .long("width")
                    .takes_value(true)
                    .min_values(1)
                    .max_values(1)
                    .multiple(true)
                    .help("Sets width of image"))
                  .arg(Arg::with_name("samples")
                    .long("samples")
                    .short("s")
                    .takes_value(true)
                    .min_values(1)
                    .max_values(1)
                    .multiple(true)
                    .help("Sets samples"))
                  .arg(Arg::with_name("threads")
                    .long("threads")
                    .short("t")
                    .takes_value(true)
                    .min_values(1)
                    .max_values(1)
                    .multiple(true)
                    .help("Set amount of worker threads"))
                  .get_matches();

  let width = if matches.is_present("width") {
    match value_t!(matches, "width", usize) {
      Ok(w) => w,
      Err(_) => 1024
    }
  } else { 1024 };

  let height = if matches.is_present("height") {
    match value_t!(matches, "height", usize) {
      Ok(h) => h,
      Err(_) => 768
    }
  } else { 768 };

  let samples = if matches.is_present("samples") {
    match value_t!(matches, "samples", usize) {
      Ok(s) => s,
      Err(_) => 2
    }
  } else { 2 };

  if matches.is_present("threads") {
    match value_t!(matches, "threads", usize) {
      Ok(t) => rayon::ThreadPoolBuilder::new().num_threads(t).build_global().unwrap(),
      Err(_) => ()
    }
  }

  let aperture = 2.6;
  let focal_length = 120.0;

  let mut events_loop = glium::glutin::EventsLoop::new();
  let window = glium::glutin::WindowBuilder::new()
    .with_dimensions((width as f64, height as f64).into())
    .with_title("rust_tracing");

  let context = glium::glutin::ContextBuilder::new()
    .with_vsync(true)
    .with_gl(GlRequest::Specific(Api::OpenGl, (3, 2)))
    .with_gl_profile(GlProfile::Core);

  let display = glium::Display::new(window, context, &events_loop).expect("Failed to create display");

  let mut buffer_texture: BufferTexture<(u8, u8, u8, u8)> =
    BufferTexture::empty_persistent(
      &display,
      width * height * 4,
      BufferTextureType::Float,
    ).expect("Failed to create rgb_buffer texture");

  {
    // init buffer texture to something
    let mut mapping = buffer_texture.map();

    for texel in mapping.iter_mut() {
      *texel = (0, 0, 0, 255);
    }
  }

  let program = glium::Program::from_source(
    &display,
    "
      #version 330 core
      void main() {
        const vec4 vertices[] = vec4[] (vec4(-1.0, -1.0, 0.5, 1.0),
                                        vec4( 1.0, -1.0, 0.5, 1.0),
                                        vec4(-1.0,  1.0, 0.5, 1.0),
                                        vec4( 1.0,  1.0, 0.5, 1.0));
        gl_Position = vertices[gl_VertexID];
      }
    ",
    "
      #version 330 core
      uniform int stride;
      uniform samplerBuffer tex;
      out vec4 color;
      void main() {
        int x = int(gl_FragCoord.x);
        int y = int(gl_FragCoord.y);
        int index = y * stride + x;
        color = texelFetch(tex, index);
      }
    ",
    None,
  ).expect("Failed to create shader");

  let mut scene_objects: Vec<Box<&dyn SceneObject>> = Vec::new();

  for tris in TRIS.iter() {
    scene_objects.push(Box::new(tris));
  }

  for sphere in SPHERES.iter() {
    scene_objects.push(Box::new(sphere));
  }

  let camera = Ray::new(Vector::new(50.0, 52.0, 295.6), Vector::new(0.0, -0.042612, -1.0).normalize());
  let focal_point = camera.org + camera.dir * focal_length;

  let cx = Vector::new(width as f64 * 0.5135 / height as f64, 0.0, 0.0);
  let cy = (cx.cross_product(&camera.dir)).normalize() * 0.5135;

  let (worker_send, main_recv) = channel::<Vec<Color>>();

  let worker = thread::spawn(move || {
    let bar = ProgressBar::new(height as u64);

    bar.set_style(ProgressStyle::default_bar()
        .template("[{elapsed_precise}] {wide_bar:.cyan/blue} {percent:>3}% {msg}"));

    // Loop over image rows.
    let total_radiances = (0..height).rev().flat_map(|y| {
      bar.inc(1);

      // Loop over image columns.
      let radiances : Vec<Color> = (0..width).into_par_iter().map(|x| {
        // 2 x 2 subsampling per pixel.
        let total_radiance = (0..2).map(|sx| {
          (0..2).map(|sy| {
            // Computes radiance at subpixel using multiple samples.
            (0..samples).map(|_| {
              // Generate random sample on circular lens.
              let random_radius = drand48();
              let random_angle = drand48();
              let lens_sample_point = aperture * Vector::new(random_radius.sqrt() * (2.0 * PI * random_angle).cos(), random_radius.sqrt() * (2.0 * PI * random_angle).sin(), 0.0);

              let dir = (focal_point - (camera.org + lens_sample_point)).normalize();

              let nu_filter_samples = || -> f64 {
                let r = 2.0 * drand48() as f64;
                if r < 1.0 { r.sqrt() - 1.0 } else { 1.0 - (2.0 - r).sqrt() }
              };

                let dx = nu_filter_samples();
                let dy = nu_filter_samples();

              // Ray direction into scene from camera through sample.
              let dir = cx * ((x as f64                       + (sx as f64 + 0.5 + dx as f64) / 2.0) / width as f64  - 0.5) +
                        cy * (((height - y - 1) as f64 + (sy as f64 + 0.5 + dy as f64) / 2.0) / height as f64 - 0.5) +
                        (camera.dir + dir).normalize();

              // Extend camera ray to start inside box.
              let start = camera.org + dir * 130.0;

              let ray = Ray::new(start + lens_sample_point, dir.normalize());

              radiance(&scene_objects, &ray, 0, 1) / samples as f64
            }).sum::<Color>().clamp(0.0, 1.0) * 0.25
          }).sum()
        }).sum();

        total_radiance
      }).collect();

      radiances
    }).collect();

    worker_send.send(total_radiances).unwrap();

    bar.finish();
  });

  let system = glium_text::TextSystem::new(&display);

  let font = glium_text::FontTexture::new(&display, &include_bytes!("../font.ttf")[..], 120, glium_text::FontTexture::ascii_character_list()).unwrap();

  let text = glium_text::TextDisplay::new(&system, &font, "Rendering...");
  let text_width = text.get_width();

  let (w, h) = display.get_framebuffer_dimensions();

  let matrix = [[2.0 / text_width, 0.0, 0.0, 0.0],
                [0.0, 2.0 * (w as f32 / 2.0) / (h as f32 / 2.0) / text_width, 0.0, 0.0],
                [0.0, 0.0, 1.0, 0.0],
                [-1.0, 0.0, 0.0, 3.0f32]];

  let mut target = display.draw();
  glium_text::draw(&text, &system, &mut target, matrix, (1.0, 1.0, 1.0, 1.0));
  target.finish().expect("failed to write text");

  let mut quit = false;

  while !quit {
    events_loop.poll_events(|event| {
      use glium::glutin::{ElementState, Event, VirtualKeyCode, WindowEvent};
      if let Event::WindowEvent { event, .. } = event {
        match event {
          WindowEvent::CloseRequested => quit = true,
          WindowEvent::KeyboardInput { input, .. } => {
            if let ElementState::Released = input.state {
              match input.virtual_keycode {
                Some(VirtualKeyCode::Escape) => quit = true,
                Some(VirtualKeyCode::F1) => {
                  let image: glium::texture::RawImage2d<u8> = display.read_front_buffer();
                  let image = image::ImageBuffer::from_raw(
                    image.width, image.height, image.data.into_owned()) .unwrap();

                  let image = image::DynamicImage::ImageRgba8(image).flipv().to_rgb();

                  image
                    .save("image.png")
                    .expect("Failed to save output image");
                },
                _ => ()
              }
            }
          }
          _ => (),
        };
      }
    });

    if quit { break; }

    if let Ok(buffer) = main_recv.recv() {
      {
        let mut mapping = buffer_texture.map();
        for (texel, rgb) in mapping.iter_mut().zip(buffer.iter()) {
          *texel = (
            (255.99 * rgb.x.min(1.0).max(0.0)) as u8,
            (255.99 * rgb.y.min(1.0).max(0.0)) as u8,
            (255.99 * rgb.z.min(1.0).max(0.0)) as u8,
            255,
          );
        }
      }

      let mut target = display.draw();

      target.draw(
        EmptyVertexAttributes { len: 4 },
        NoIndices(PrimitiveType::TriangleStrip),
        &program,
        &uniform!{ tex: &buffer_texture, stride: width as i32 },
        &Default::default(),
      ).unwrap();

      target.finish().unwrap();
    }
  }

  worker.join().unwrap();
}
