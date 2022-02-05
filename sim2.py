from math import log10
from classes import proton, electron
from vpython import vector, rate

dt = 1e-23


def update_forces(particles):
    closest = 1e99
    for part in particles:
        part.set_force(particles)
        if part.closest < closest:
            closest = part.closest
    return closest


def update(particles):
    global dt
    closest = update_forces(particles)
    for part in particles:
        part.set_velocity(dt)
        part.set_pos(dt)
        print(particles.index(part), " : ", part.force,
              " : ", part.velocity, " : ", part.pos)

    print(dt)
    if closest > 0:
        dt = 10**(-7.15743510167 * 1.08695652174**(-log10(closest)))


def main():
    particles = []
    for i in range(2):
        particles.append(proton(identity=i, pos=vector(
            i*-5e-14, 0, 0), velocity=vector(2500000 * (-1)**(i+1), 0, 0)))
    particles.append(electron(identity=2, pos=vector(
        -5e-14, 1e-12, 0), velocity=vector(500000, 50000, 200000)))

    while True:
        rate(200)
        # input()
        update(particles)


if __name__ == "__main__":
    main()
