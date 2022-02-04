from vpython import vector, rate, sqrt, scene, pi, sphere
from math import hypot
# GlowScript 3.2 VPython

SCALE = 1e16
G = SCALE * 6.67384e-11  # Nm^2kg^-2
f = SCALE * 8.987551787e9  # Nm^2C^-2
m_proton = SCALE * 1.672621777e-27  # kg
m_neutron = SCALE * 1.674927351e-27  # kg
m_electron = SCALE * 9.10938291e-31  # kg
e = SCALE * 1.602176565e-19  # C

radius_proton = radius_neutron = SCALE * 8.4e-16  # m
radius_electron = radius_proton / 2  # m
radius_sim = radius_proton * 1e15 / 5
balls = [sphere(make_trail=True, trail_type="points", radius=radius_sim), sphere(
    make_trail=True, trail_type="points", radius=radius_sim)]
for index, ball in enumerate(balls):
    ball.pos = vector(index * 1e15 * radius_proton, 0, 0)
    ball.velocity = vector(0, 0, 0)
    ball.force = vector(0, 0, 0)
    ball.mass = m_proton
    ball.charge = e


def update_particle(b, t) -> None:
    b.force = vector(0, 0, 0)
    self_ind = balls.index(b)
    for i, ball_elem in enumerate(balls):
        if i != self_ind:
            r_vec = vector(b.pos.x-ball_elem.pos.x, b.pos.y -
                           ball_elem.pos.y, b.pos.z-ball_elem.pos.z)
            r = hypot(b.pos.x-ball_elem.pos.x, b.pos.y -
                      ball_elem.pos.y, b.pos.z-ball_elem.pos.z)
            r2 = r**2
            print(r2)
            f_g = G*b.mass * ball_elem.mass / r2
            b.force += r_vec * f_g / r
            f_em = f*b.charge * ball_elem.charge / r2
            b.force += r_vec * f_em / r


def update(t):
    for ball_elem in balls:
        update_particle(ball_elem, t)
    for ball_elem2 in balls:
        ball_elem2.velocity += ball_elem2.force / ball_elem2.mass
        print(balls.index(ball_elem2), ": ",
              ball_elem2.force, ": ", ball_elem2.velocity)
        ball_elem2.pos += ball_elem2.velocity


dt = 0.1
time = 0
while True:
    rate(1)
    update(time)
    time += dt
