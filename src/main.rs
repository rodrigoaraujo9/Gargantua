use physical_constants::{NEWTONIAN_CONSTANT_OF_GRAVITATION, SPEED_OF_LIGHT_IN_VACUUM};
use raylib::prelude::*;

pub const SCREEN_WIDTH: i32 = 800;
pub const SCREEN_HEIGHT: i32 = 600;

pub const G: f64 = NEWTONIAN_CONSTANT_OF_GRAVITATION;
pub const C: f64 = SPEED_OF_LIGHT_IN_VACUUM; //m/s
pub const SOLAR_MASS: f64 = 1.988416e30f64; //kg
pub const SCALE_FACTOR: f64 = 3.4e-11f64; //scale down
pub const SPEEDUP_FACTOR: f64 = 39235.0; // acceleration of frame so light movement is visible -> if 1, real time

struct Beem {
    position: Vector2,
    angle: f32, // 0 - 2PI
}

impl Beem {
    const RADIUS: f32 = 5.0;
    pub fn new(x: f32, y: f32, angle: f32) -> Self {
        Self {
            position: Vector2::new(x, y),
            angle,
        }
    }

    fn update(&mut self, d: &mut RaylibDrawHandle) {
        let distance: f32 = (d.get_frame_time() as f64 * C * SCALE_FACTOR * SPEEDUP_FACTOR) as f32;
        self.position.x += distance * self.angle.cos();
        self.position.y += distance * self.angle.sin();

        if self.position.x > SCREEN_WIDTH as f32 {
            self.position.x = 0.0;
        }
        if self.position.x < 0.0 {
            self.position.x = SCREEN_WIDTH as f32;
        }
        if self.position.y > SCREEN_HEIGHT as f32 {
            self.position.y = 0.0;
        }
        if self.position.y < 0.0 {
            self.position.y = SCREEN_HEIGHT as f32;
        }
    }

    //collision detection
    fn within_event_horizon(&self, black_hole: &BlackHole) -> bool {
        let distance = ((self.position.x - black_hole.position.x).powi(2)
            + (self.position.y - black_hole.position.y).powi(2))
        .sqrt();
        distance <= black_hole.event_horizon_radius() as f32
    }

    pub fn draw(&self, d: &mut RaylibDrawHandle) {
        d.draw_circle_v(self.position, Self::RADIUS, Color::WHITE);
    }
}

impl Drop for Beem {
    fn drop(&mut self) {
        // destructor
    }
}

struct BlackHole {
    position: Vector2,
    mass: f64,   //in solar masses
    radius: f64, // in pixels
}

impl BlackHole {
    pub fn new(x: f32, y: f32, mass: f64) -> Self {
        Self {
            position: Vector2::new(x, y),
            mass,
            radius: Self::event_horizon_radius_init(mass),
        }
    }
    // public r_s
    pub fn event_horizon_radius(&self) -> f64 {
        let mass_kg: f64 = self.mass * SOLAR_MASS;
        let r_s: f64 = SCALE_FACTOR * (2.0 * G * mass_kg) / (C * C);
        r_s
    }

    // r_s used for init
    fn event_horizon_radius_init(mass: f64) -> f64 {
        let mass_kg: f64 = mass * SOLAR_MASS;
        let r_s: f64 = SCALE_FACTOR * (2.0 * G * mass_kg) / (C * C);
        r_s
    }

    pub fn draw(&self, d: &mut RaylibDrawHandle) {
        let r_s = self.event_horizon_radius();
        d.draw_circle_v(self.position, r_s as f32, Color::BLACK);
    }
}

fn main() {
    let (mut rl, thread) = raylib::init()
        .size(SCREEN_WIDTH, SCREEN_HEIGHT)
        .title("Gargantula")
        .build();
    rl.set_target_fps(60);
    let sh: f32 = SCREEN_HEIGHT as f32;
    let sw: f32 = SCREEN_WIDTH as f32;
    let gargantula = BlackHole::new(sw / 2.0, sh / 2.0, 10e8);
    let mut beem: Option<Beem> = Some(Beem::new(0.0, sh / 2.0, 0.0));

    while !rl.window_should_close() {
        let mut d = rl.begin_drawing(&thread);
        d.clear_background(Color::LIGHTGRAY);
        if let Some(ref mut beam) = beem {
            beam.update(&mut d);
            beam.draw(&mut d);
            if beam.within_event_horizon(&gargantula) {
                beem = None;
            }
        }
        gargantula.draw(&mut d);

        let mass_text = format!("mass: {:.1e} solar masses", gargantula.mass);
        let radius_km = gargantula.event_horizon_radius() / SCALE_FACTOR;
        let radius_text = format!("event horizon: {:.1} km", radius_km);
        let speedup_text = format!("speedup: x{:.1e}", SPEEDUP_FACTOR);

        d.draw_text(&mass_text, 10, 10, 20, Color::BLACK);
        d.draw_text(&radius_text, 10, 35, 20, Color::BLACK);
        d.draw_text(&speedup_text, 10, 60, 20, Color::BLACK);
    }
}
