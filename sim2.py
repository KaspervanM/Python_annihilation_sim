from classes import proton
from vpython import vector, rate, hypot


def update_forces(particles):
    for part in particles:
        part.set_force(particles)


def update(particles, dt):
    update_forces(particles)
    for part in particles:
        part.set_velocity()
        part.set_pos(dt)
        print(particles.index(part), " : ",
              part.force, " : ", part.velocity, " : ", part.pos)


def main():
    particles = []
    for i in range(2):
        particles.append(proton(identity=i, pos=vector(
            i*1e-1, 0, 0), velocity=vector(50 * (-1)**i, 0, 0)))

    dt = 0.00001
    while True:
        rate(1)
        update(particles, dt)


if __name__ == "__main__":
    main()
