import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import RegularGridInterpolator
from scipy import interpolate
import math,sys,os,pathlib

nbnds=120
nkx,nky=71,71
kxc,kyc=0.,0.
kxmax,kymax=0.05,0.05
ef=-1.1234
ef_shift=0.87
fn='sbp.dat'
spin_direction='z'
plt_2Dfermi=True
plt_arrow=False
plt_spol=True

plt.rcParams['font.family'] = 'Times New Roman'
plt.rcParams["xtick.labelsize"]= 15.0
plt.rcParams["ytick.labelsize"]= 15.0
plt.rcParams["xtick.major.pad"] = 5
plt.rcParams["ytick.major.pad"] = 5
plt.rcParams["axes.labelsize"] = 15.0
plt.rcParams["axes.linewidth"] = 1.0
plt.rcParams["axes.labelpad"] = 6
plt.rcParams["xtick.direction"] = "in" 
plt.rcParams["ytick.direction"] = "in"
plt.rcParams["xtick.major.width"] = 1.0
plt.rcParams["ytick.major.width"] = 1.0
plt.rcParams["xtick.minor.width"] = 0.5
plt.rcParams["ytick.minor.width"] = 0.5
plt.rcParams["xtick.major.size"] = 4.5
plt.rcParams["ytick.major.size"] = 4.5
plt.rcParams["xtick.minor.size"] = 3.0
plt.rcParams["ytick.minor.size"] = 3.0

def mk_kgrid(nkx,nky,kxc=0.,kyc=0.,kxmax=0.02,kymax=0.02,kz=0.,code='qe'):
    kx=np.linspace(kxc-kxmax,kxc+kxmax,nkx)
    ky=np.linspace(kyc-kymax,kyc+kymax,nky)
    if code=='qe':
        kgrids=f"{nkx*nky} \n"+"\n".join([f" {kxx}  {kyy}  {kz}  1 " for kxx in kx for kyy in ky])
        with open("Kpoints-2D.txt","w") as f:
            f.write(kgrids)
    print("* 2D k-points are written out Kpoints-2D.txt!")
    return kx,ky

def read_eigs(fn,nbnds,nkx,nky,code='qe'):
    nks=nkx*nky
    if code=='qe':
        nl=math.ceil(float(nbnds)/10.)
        nr=nbnds%10
        eig=[]
        if nr==0:
            nr=10
        ic=1
        with open(fn,'r') as f:
            data=f.readlines()
        for dat in data:
            dat=dat.split()
            if ic==1:
                pass
            elif (ic-2)%(nl+1)==0:
                pass
            elif dat[0]=='\n':
                pass
            elif (ic-1)%(nl+1)==0:
                for i in range(nr):
                    eig.append(float(dat[i]))
            else:
                for i in range(10):
                    eig.append(float(dat[i]))
            ic+=1
        eig=np.array(eig).reshape(nkx,nky,nbnds)
    print(f"* read eigenvalues in {fn}")
    print(eig)
    return eig
                
def plt_2Dfs(kx,ky,eig,sx,sy,sz,nbnds,ef=0.,efsh=0.,kz=0.,kcx=0.,kcy=0.,kxmax=0.02,kymax=0.02,sd='z',
             plt_2Dfermi=False,plt_arrow=False,plt_spol=False,tol=1e-3,nd=501):
    plt.figure(figsize=(8,6))
    plt.xlabel(r'$k_x$')
    plt.ylabel(r'$k_y$')
    plotted = False  # Flag to track if anything was plotted
    for i in range(nbnds):
        if eig[:, :, i].max() < ef + efsh or eig[:, :, i].min() > ef + efsh:
            continue
        eigip = RegularGridInterpolator((kx, ky), eig[:, :, i])
        sxip = RegularGridInterpolator((kx, ky), sx[:, :, i])
        syip = RegularGridInterpolator((kx, ky), sy[:, :, i])
        szip = RegularGridInterpolator((kx, ky), sz[:, :, i])
        if sd == 'x':
            spinip = sxip
            spinip1 = syip
            spinip2 = szip
        elif sd == 'y':
            spinip = syip
            spinip1 = szip
            spinip2 = sxip
        else:  # 'z'
            spinip = szip
            spinip1 = sxip
            spinip2 = syip
        kxd,kyd=mk_kgrid(nd,nd,kcx=kcx,kcy=kcy,kxmax=kxmax,kymax=kymax)
        kxd,kyd = np.meshgrid(kxd,kyd,indexing='ij')
        if plt_2Dfermi:
            cs = plt.contour(kxd,kyd,eigip((kxd,kyd)),levels=[ef + efsh],
                             colors=['lightgray'], linestyles=['solid'])
            if not cs.collections or len(cs.collections[0].get_paths()) == 0:
                continue  # Skip if no contour found
            for path in cs.collections[0].get_paths():
                v = path.vertices
                plotted = True
                if plt_spol:
                    valid_mask = (
                        (v[:, 0] >= kx.min()) & (v[:, 0] <= kx.max()) &
                        (v[:, 1] >= ky.min()) & (v[:, 1] <= ky.max())
                    )
                    v_valid = v[valid_mask]
                    if len(v_valid) > 0:
                        fig = plt.scatter(v_valid[:, 0], v_valid[:, 1], vmin=-0.5, vmax=0.5,
                                          marker="o", c=spinip(v_valid), cmap='bwr')
                if plt_arrow:
                    nskip = 25
                    plt.quiver(v[::nskip, 0], v[::nskip, 1], spinip1(v[::nskip, :]),
                               spinip2(v[::nskip, :]), color='black', angles='xy',
                               scale_units='xy', scale=100, width=5e-3)
    if plt_spol and plotted and 'fig' in locals():
        cbar = plt.colorbar(fig, label=f'$s_{sd}$')
        cbar.set_ticks([-0.5, 0., 0.5])
    if plotted and filename:
        plt.savefig(filename, dpi=300)
        print(f"* Saved plot to {filename}")
    else:
        print(f"- No Fermi surface for ef_shift = {efsh:.2f}, skipping save.")
    plt.close()
    
if __name__=='__main__':
    kx,ky=mk_kgrid(nkx,nky,kcx=kcx,kcy=kcy,kxmax=kxmax,kymax=kymax)
    if not (plt_2Dfermi or plt_arrow or  plt_spol):
        sys.exit()
    eig=read_eigs(fn,nbnds,nkx,nky)
    sx=read_eigs(fn+".1",nbnds,nkx,nky)
    sy=read_eigs(fn+".2",nbnds,nkx,nky)
    sz=read_eigs(fn+".3",nbnds,nkx,nky)
    plt_2Dfs(kx,ky,eig,sx,sy,sz,nbnds,ef=ef,efsh=ef_shift,kz=0.,kcx=kcx,kcy=kcy,kxmax=kxmax,kymax=kymax,sd=spin_direction,
             plt_2Dfermi=plt_2Dfermi,plt_arrow=plt_arrow,plt_spol=plt_spol,tol=3e-2)
