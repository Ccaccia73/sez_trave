import numpy as np

class Section:
    def __init__(self, a, b, shape):
        self.a = a
        self.b = b
        self.shape = shape
        
    
    def setup(self, na, nb, grade, E, nu, off_x=0., off_y=0.):
        self.na = na
        self.nb = nb
        self.gr = grade
        self.E = E
        self.nu = nu
        self.G = self.E/(2*(1+self.nu))
        self.ofsx = off_x
        self.ofsy = off_y
    
    def generazione_dati(self):
        if self.gr == 2:
            self.xx = np.linspace(-1., 1., num=self.na+1, endpoint=True)
            self.yy = np.linspace(-1., 1., num=self.nb+1, endpoint=True)
            
            self.ib = np.array([1., 2., 1+(self.nb*(self.gr-1)+1), 2+(self.nb*(self.gr-1)+1)])
            self.ad = np.ones(4)
            
        elif self.gr == 3:
            self.xx = np.linspace(-1., 1., num=2*self.na+1, endpoint=True)
            self.yy = np.linspace(-1., 1., num=2*self.nb+1, endpoint=True)
            
            self.ib = np.array([1., 2., 3., \
            1+(self.nb*(self.gr-1)+1), 2+(self.nb*(self.gr-1)+1), 3+(self.nb*(self.gr-1)+1), \
            1+(self.nb*(self.gr-1)+1)*2, 2+(self.nb*(self.gr-1)+1)*2, 3+(self.nb*(self.gr-1)+1)*2])
            
            self.ad = 2*np.ones(9)
            
        self.tabe = np.zeros((self.na*self.nb, self.gr**2))
        
        ie = 0
        num0 = self.ib
        
        for iel in range(self.na):
            num = num0+iel*(self.nb*(self.gr-1)+1)*(self.gr-1)
            
            for jel in range(self.nb):
                self.tabe[ie,:] = num
                ie+=1
                num+=self.ad
        
        nx = len(self.xx)
        ny = len(self.yy)

        self.nt = nx*ny

        self.xnd = np.tile(self.xx,(self.na*(self.gr-1)+1,1)).reshape((self.na*(self.gr-1)+1)**2 ,order='F')
        self.ynd = np.tile(self.yy,self.na*(self.gr-1)+1)
        
        self.F = np.zeros((6,6))
        
        self.F[0,0] = 1./self.E
        self.F[1,1] = self.F[0,0]
        self.F[5,5] = self.F[0,0]
        self.F[2,2] = 1/self.G
        self.F[3,3] = self.F[2,2]
        self.F[4,4] = self.F[2,2]
        
        self.F[0,1] = self.nu/self.E
        self.F[0,2] = self.F[0,1]
        self.F[0,5] = self.F[0,1]
        self.F[1,0] = self.F[0,1]
        self.F[1,5] = self.F[0,1] 
        self.F[5,0] = self.F[0,1]
        self.F[5,1] = self.F[0,1]
        
        self.D = np.linalg.inv(self.F)
        t = [2, 0, 1, 3, 4, 5]
        # http://docs.scipy.org/doc/numpy/reference/arrays.indexing.html#advanced-indexing
        self.D = self.D[:, t][t]
        
        if self.shape == 'rect':
            self.area = self.a*self.b
            self.I1 = 1/12*self.b**3*self.a
            self.I2 = 1/12*self.a**3*self.b
            
            ratio = self.a / self.b            
            
            if ratio >= 1:
                self.Jt = self.a*self.b**3/16*(16/3-3.61558*ratio*(1-ratio**4/12))
            else:
                self.Jt = self.b*self.a**3/16*(16/3-3.61558*(1/ratio)*(1-(1/ratio)**4/12))
                        
            self.xnd = self.xnd*self.a/2 + self.ofsx
            self.ynd = self.ynd*self.b/2 + self.ofsy
            
        elif self.shape == 'ell' or self.shape == 'cir':
            
            a1 = self.a
            if self.shape == 'ell':
                b1 = self.b
            else:
                b1 = a1
            
            self.area = np.pi*a1*b1
            self.I1=np.pi*b1**3*a1/4
            self.I2=np.pi*a1**3*b1/4
            self.Jt=np.pi*a1**3*b1**3/(a1**2+b1**2)
            
            cs=np.cos(45*np.pi/180)
            XG=np.array([-cs, -1, -cs, 0, 0, 0, cs, 1, cs])*a1;
            YG=np.array([-cs,  0, cs, -1,  0, 1, -cs,  0, cs])*b1;

            for i_n, (r,s) in enumerate(zip(self.xnd, self.ynd)):
                
                H1dr=np.array([1/2*(1-r)-1/2*(1-r**2), 1-r**2, 1/2*(1+r)-1/2*(1-r**2)])
                H1ds=np.array([1/2*(1-s)-1/2*(1-s**2), 1-s**2, 1/2*(1+s)-1/2*(1-s**2)])
                
                H = np.outer(H1ds, H1dr)
                
                xnn = np.dot(H,XG)
                ynn = np.dot(H,YG)
                
                
                if abs((xnn/a1)**2+(ynn/b1)**2 -1) < 0.06:
                    if self.shape == 'cir':
                        ang=np.arctan2(ynn,xnn)
                        xnn=np.cos(ang)*a1
                        ynn=np.sin(ang)*b1
                    else:
                        ynn = np.sign(ynn) * np.sqrt((1-(xnn/a1)**2) *b1**2 )
                            
                    
                self.xnd[i_n]=xnn+self.ofsx
                self.ynd[i_n]=ynn+self.ofsy
                
        else:
            print("Forma sconosciuta")

        self.kroark = np.array([self.G*self.area*5/6, self.G*self.area*5/6, self.E*self.area, self.E*self.I1, self.E*self.I2, self.G*self.Jt])
        
        
    def assemblaggio(self):
        
        self.AG=np.zeros((6,6))
        self.RG=np.zeros((6,self.nt*3))
        self.LG=np.zeros((6,self.nt*3))
        self.EG=np.zeros((self.nt*3,self.nt*3))
        self.CG=np.zeros((self.nt*3,self.nt*3))
        self.MG=np.zeros((self.nt*3,self.nt*3))
        
        
        
        
