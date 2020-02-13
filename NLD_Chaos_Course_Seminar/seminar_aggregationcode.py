import pygame,random, math
from pygame.locals import *

xmax = 900    #width of window
ymax = 900     #height of window
v_limit= 500
radius = 4
maxclst = 1
n_particles= 800

class Particle:

   def __init__(self, startx, starty, col):
       self.x = startx
       self.y = starty
       self.col = col
       self.sx = startx
       self.sy = starty
       self.clst_id=0
       self.clst_no=1



   def move(self,angle):
       #a = random.uniform(0, 2*math.pi)
       if self.clst_id == 0:
           angle= random.uniform(0, 2*math.pi)
       self.y+= math.sin(angle)* v_limit/self.clst_no
       self.x+= math.cos(angle)* v_limit/self.clst_no
       if self.y < 0:
           #self.x = self.sx
           self.y += 900
       if self.y > 900:
           self.y-=900
       if self.x < 0:
               # self.x = self.sx
               self.x+= 900
       if self.x > 900:
               self.x -= 900
           #self.y -= 1

   def check_stick(self, Particle):
       global maxclst
       dist = math.hypot(self.x-Particle.x, self.y - Particle.y)
       if dist < 2*radius : # to ckeck if it sticks
           if self.clst_id == 0  and  Particle.clst_id ==0 :
             self.clst_id = maxclst
             Particle.clst_id = maxclst
             maxclst += 1

           elif self.clst_id > Particle.clst_id :
               Particle.clst_id = self.clst_id

           else :
               self.clst_id = Particle.clst_id

           self.clst_no += 1
           Particle.clst_no += 1


def main():

   pygame.init()
   screen = pygame.display.set_mode((xmax,ymax))
   white = (255, 255, 255)
   black = (0,0,0)
   #grey = (128,128,128)

   clock=pygame.time.Clock()

   particles = []
   cluster_angles= []
   for part in range(n_particles):
       col = black
       particles.append(Particle(random.randint(0,900), random.randint(0,900), col))

   exitflag = False
   while not exitflag:
       for event in pygame.event.get():
           if event.type == QUIT:
               exitflag = True
           elif event.type == KEYDOWN:
               if event.key == K_ESCAPE:
                   exitflag = True

       cluster_angles[:]= []
       for i in range(0,n_particles):
           cluster_angles.append(random.uniform(0, 2*math.pi))

       screen.fill(white)
       for p in particles:
           p.move(cluster_angles[p.clst_id])
           p.x = int(p.x)
           p.y = int(p.y)
           pygame.draw.circle(screen, p.col, (p.x, p.y), radius)
       #for p1 in particles:
           for p2 in particles:
          	 if p2!=p: 	
           		p.check_stick(p2)
           		
       pygame.display.flip()
       clock.tick(500)
   pygame.quit()

if __name__ == "__main__":
   main()

