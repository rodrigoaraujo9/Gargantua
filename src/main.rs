use physical_constants::{NEWTONIAN_CONSTANT_OF_GRAVITATION, SPEED_OF_LIGHT_IN_VACUUM};
use raylib::prelude::*;
use std::f64::consts::PI;

pub const SCREEN_WIDTH: i32 = 800;
pub const SCREEN_HEIGHT: i32 = 600;
pub const NUM_BEEMS: i32 = 1000;

pub const G: f64 = NEWTONIAN_CONSTANT_OF_GRAVITATION;
pub const C: f64 = SPEED_OF_LIGHT_IN_VACUUM; //m/s
pub const SOLAR_MASS: f64 = 1.988416e30f64; //kg
pub const SCALE_FACTOR: f64 = 2.0 * 3.4e-11f64; //scale down
pub const SPEEDUP_FACTOR: f64 = 0.0000015; // acceleration of frame so light movement is visible

// angular step size (dphi) for simulation speed
pub const PHI_STEP_FACTOR: f64 = 0.01;

#[derive(Copy, Clone)]
struct BeemState {
    u: f64,
    phi: f64,
    du_dphi: f64,
}

struct Beem {
    position: Vector2, // position in screen coordinates
    angle: f32,        // 0 - 2PI
    trail: Vec<Vector2>,
    is_alive: bool,

    // relevant for respawn logic
    original_position: Vector2,
    original_angle: f32,

    u: f64,       // 1/r -> inverse of the physical distance from the black hole
    phi: f64,     // angular position
    du_dphi: f64, // angular momentum

    // flag ensure physics are setup only once
    initialized: bool,
}

impl Beem {
    const RADIUS: f32 = 2.0;
    const MAX_TRAIL_LENGTH: usize = 100;

    pub fn new(x: f32, y: f32, angle: f32) -> Self {
        Self {
            position: Vector2::new(x, y),
            angle,
            trail: Vec::new(),
            is_alive: true,
            original_position: Vector2::new(x, y),
            original_angle: angle,

            u: 0.0,
            phi: 0.0,
            du_dphi: 0.0,
            initialized: false,
        }
    }

    fn initialize_geodesics(&mut self, black_hole: &BlackHole) {
        // once per beem
        if self.initialized {
            return;
        }

        // convert screen coordinates to physical polar coordinates relative to the black hole
        let dx = self.position.x - black_hole.position.x;
        let dy = self.position.y - black_hole.position.y;
        let r = ((dx * dx + dy * dy) as f64).sqrt() / SCALE_FACTOR;
        let phi = (dy as f64).atan2(dx as f64);

        self.u = 1.0 / r;
        self.phi = phi;

        // calculate the initial angular momentum from the beems initial direction
        let angle_to_bh = (dy as f64).atan2(dx as f64);
        let relative_angle = self.angle as f64 - angle_to_bh;

        self.du_dphi = -relative_angle.tan() * (1.0 / r);

        // flag as init to prevent redoing setup
        self.initialized = true;
    }

    fn get_derivatives(&self, state: BeemState, black_hole: &BlackHole) -> (f64, f64) {
        // calulate the derivatives from the geodesic equation
        let r_s_physical = black_hole.event_horizon_radius_physical();
        let d2u_dphi2 = -state.u + (1.5 * r_s_physical) * state.u * state.u;
        (state.du_dphi, d2u_dphi2)
    }

    // rk4 method for calculating physics, averaging 4 steps to predict physics more acurately
    fn rk4_step(&self, state: BeemState, dphi: f64, black_hole: &BlackHole) -> BeemState {
        let (k1_u, k1_du_dphi) = self.get_derivatives(state, black_hole);

        let k2_state = BeemState {
            u: state.u + k1_u * dphi / 2.0,
            phi: state.phi + dphi / 2.0,
            du_dphi: state.du_dphi + k1_du_dphi * dphi / 2.0,
        };
        let (k2_u, k2_du_dphi) = self.get_derivatives(k2_state, black_hole);

        let k3_state = BeemState {
            u: state.u + k2_u * dphi / 2.0,
            phi: state.phi + dphi / 2.0,
            du_dphi: state.du_dphi + k2_du_dphi * dphi / 2.0,
        };
        let (k3_u, k3_du_dphi) = self.get_derivatives(k3_state, black_hole);

        let k4_state = BeemState {
            u: state.u + k3_u * dphi,
            phi: state.phi + dphi,
            du_dphi: state.du_dphi + k3_du_dphi * dphi,
        };
        let (k4_u, k4_du_dphi) = self.get_derivatives(k4_state, black_hole);

        BeemState {
            u: state.u + (k1_u + 2.0 * k2_u + 2.0 * k3_u + k4_u) * dphi / 6.0,
            phi: state.phi + dphi,
            du_dphi: state.du_dphi
                + (k1_du_dphi + 2.0 * k2_du_dphi + 2.0 * k3_du_dphi + k4_du_dphi) * dphi / 6.0,
        }
    }

    fn step(&mut self, black_hole: &BlackHole, r_s: f64, dphi: f64) {
        // make use of rk4
        let current_state = BeemState {
            u: self.u,
            phi: self.phi,
            du_dphi: self.du_dphi,
        };

        let new_state = self.rk4_step(current_state, dphi, black_hole);

        // update the light beem with new phys state
        self.u = new_state.u;
        self.phi = new_state.phi;
        self.du_dphi = new_state.du_dphi;

        let r = 1.0 / self.u;

        // check if array has fallen beyond the event horizon
        if r < r_s * 1.01 {
            self.is_alive = false;
            return;
        }

        // convert back to screen coordinates for drawing -> cartesian
        let x = (r * self.phi.cos()) * SCALE_FACTOR + black_hole.position.x as f64;
        let y = (r * self.phi.sin()) * SCALE_FACTOR + black_hole.position.y as f64;

        self.position = Vector2::new(x as f32, y as f32);
    }

    // handles state for beem
    fn update(&mut self, d: &mut RaylibDrawHandle, black_hole: &BlackHole) {
        // fades the trail if dead
        if !self.is_alive {
            if !self.trail.is_empty() {
                self.trail.remove(0);
            }
            return;
        }

        // initializes the beems physical state once and only once
        if !self.initialized {
            self.initialize_geodesics(black_hole);
        }

        self.trail.push(self.position);

        // limit trail lenght for fading effect
        if self.trail.len() > Self::MAX_TRAIL_LENGTH {
            self.trail.remove(0);
        }

        // check if beem is beyond event horizon
        let r_s_physical = black_hole.event_horizon_radius_physical();
        let r = 1.0 / self.u;

        if r < r_s_physical {
            self.is_alive = false;
            return;
        }

        // advance the beem's state using the geodesic equations
        let dphi = PHI_STEP_FACTOR; // Use the new constant
        self.step(black_hole, r_s_physical, dphi);

        // kills the ray if it goes off-screen
        if self.position.x > SCREEN_WIDTH as f32
            || self.position.x < 0.0
            || self.position.y > SCREEN_HEIGHT as f32
            || self.position.y < 0.0
        {
            self.is_alive = false;
        }
    }

    // auxilary functions
    pub fn should_remove(&self) -> bool {
        !self.is_alive && self.trail.is_empty()
    }

    pub fn respawn(&mut self) {
        self.position = self.original_position;
        self.angle = self.original_angle;
        self.trail.clear();
        self.is_alive = true;

        self.u = 0.0;
        self.phi = 0.0;
        self.du_dphi = 0.0;
        self.initialized = false;
    }

    // draw funct that handles the liberties taken to represent the light beem
    pub fn draw(&self, d: &mut RaylibDrawHandle) {
        for i in 1..self.trail.len() {
            let alpha = (i as f32 / self.trail.len() as f32 * 255.0) as u8;
            let trail_color = Color::new(169, 169, 169, alpha);
            d.draw_line_ex(
                self.trail[i - 1],
                self.trail[i],
                Self::RADIUS * 2.0,
                trail_color,
            );
        }

        if self.is_alive {
            d.draw_circle_v(self.position, Self::RADIUS, Color::WHITE);
        }
    }
}

// destructor logic for beem
impl Drop for Beem {
    fn drop(&mut self) {}
}

struct BlackHole {
    position: Vector2,
    mass: f64,
    radius: f64,
}

impl BlackHole {
    pub fn new(x: f32, y: f32, mass: f64) -> Self {
        Self {
            position: Vector2::new(x, y),
            mass,
            radius: Self::event_horizon_radius_init(mass),
        }
    }

    pub fn event_horizon_radius_physical(&self) -> f64 {
        (2.0 * G * (self.mass * SOLAR_MASS)) / (C * C)
    }

    pub fn event_horizon_radius(&self) -> f64 {
        SCALE_FACTOR * self.event_horizon_radius_physical()
    }

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
        .title("Gargantua")
        .build();
    rl.set_target_fps(60);
    let sh: f32 = SCREEN_HEIGHT as f32;
    let sw: f32 = SCREEN_WIDTH as f32;

    let gargantua = BlackHole::new(sw / 2.0, sh / 2.0, 10e8);

    let mut beems: Vec<Beem> = Vec::new();
    let num_beams = NUM_BEEMS / 4;
    let beam_spacing = sh / (num_beams + 1) as f32;

    for i in 1..=num_beams {
        let y_pos = i as f32 * beam_spacing;
        beems.push(Beem::new(0.0, y_pos, 0.0));
    }

    for i in 1..=num_beams {
        let y_pos = i as f32 * beam_spacing;
        let lim: f32 = (SCREEN_WIDTH - 1) as f32;
        beems.push(Beem::new(lim, y_pos, PI as f32));
    }

    let beam_spacing_x = sw / (num_beams + 1) as f32;

    for i in 1..=num_beams {
        let x_pos = i as f32 * beam_spacing_x;
        beems.push(Beem::new(x_pos, 0.0, (PI as f32) / 2.0));
    }

    for i in 1..=num_beams {
        let x_pos = i as f32 * beam_spacing_x;
        let lim2: f32 = (SCREEN_HEIGHT - 1) as f32;
        beems.push(Beem::new(x_pos, lim2, 3.0 * (PI as f32) / 2.0));
    }

    while !rl.window_should_close() {
        let mut d = rl.begin_drawing(&thread);

        d.clear_background(Color::new(20, 20, 20, 255));

        for beem in &mut beems {
            beem.update(&mut d, &gargantua);
            beem.draw(&mut d);

            if beem.should_remove() {
                beem.respawn();
            }
        }

        gargantua.draw(&mut d);

        let mass_text = format!("mass: {:.1e} solar masses", gargantua.mass);
        let light_beems_text = format!("number of light beems: {:}", NUM_BEEMS);
        let radius_km = gargantua.event_horizon_radius_physical() / 1000.0;
        let radius_text = format!("event horizon: {:.1} km", radius_km);
        let speedup_text = format!("angular step: {}", PHI_STEP_FACTOR);

        d.draw_text(&mass_text, 10, 10, 20, Color::WHITE);
        d.draw_text(&radius_text, 10, 35, 20, Color::WHITE);
        d.draw_text(&light_beems_text, 10, 60, 20, Color::WHITE);
        d.draw_text(&speedup_text, 10, 85, 20, Color::WHITE);
    }
}
