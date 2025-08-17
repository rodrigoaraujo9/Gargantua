use physical_constants::{NEWTONIAN_CONSTANT_OF_GRAVITATION, SPEED_OF_LIGHT_IN_VACUUM};
use raylib::prelude::*;

pub const SCREEN_WIDTH: i32 = 800;
pub const SCREEN_HEIGHT: i32 = 600;

struct BlackHole {
    position: Vector2,
    mass: f64, //in solar masses
}
impl BlackHole {
    pub const G: f64 = NEWTONIAN_CONSTANT_OF_GRAVITATION;
    pub const C: f64 = SPEED_OF_LIGHT_IN_VACUUM; //m/s
    pub const SOLAR_MASS: f64 = 1.988416e30f64; //kg
    pub const SCALE_FACTOR: f64 = 3.4e-11f64; //scale down
    pub fn new(x: f32, y: f32, mass_: f64) -> Self {
        Self {
            position: Vector2::new(x, y),
            mass: mass_,
        }
    }
    // r_s
    fn event_horizon_radius(&self) -> f64 {
        let mass_kg: f64 = self.mass * Self::SOLAR_MASS;
        let r_s: f64 = Self::SCALE_FACTOR * (2.0 * Self::G * mass_kg) / (Self::C * Self::C);
        r_s
    }

    fn draw(&self, d: &mut RaylibDrawHandle) {
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

    while !rl.window_should_close() {
        let mut d = rl.begin_drawing(&thread);
        d.clear_background(Color::WHITE);
        gargantula.draw(&mut d);
    }
}
