function MOut=BlanchMat(MIn,F,dt)

%MOut=BlanchMat(MIn,F,dt)

TF = fft(MIn);

npts = size(MIn,1);
dw   = 1/(npts*dt);

Nmin = floor(F(1)/dw)+1;
Nmax = floor(F(2)/dw)+1;

Napod=min([20 floor((Nmax-Nmin+1)/2)]);

N1=max([1 Nmin-Napod]);
N2=Nmin+Napod;
N3=Nmax-Napod;
N4=min([floor((npts-1)/2) Nmax+Napod]);

N5=max([npts - Nmax - Napod floor((npts-1)/2)]);
N6=npts - Nmax + Napod;
N7=npts - Nmin - Napod;
N8=min([npts - Nmin + Napod npts]);

MOut=zeros(size(MIn));

MOut(N1:N2,:)=repmat(cos((pi/2:-pi/(2*(N2-N1)):0)).',1,size(MIn,2));
MOut(N2:N3,:)=1;
MOut(N3:N4,:)=repmat(cos((0:pi/(2*(N4-N3)):pi/2)).',1,size(MIn,2));

MOut(N5:N6,:)=repmat(cos((pi/2:-pi/(2*(N6-N5)):0)).',1,size(MIn,2));
MOut(N6:N7,:)=1;
MOut(N7:N8,:)=repmat(cos((0:pi/(2*(N8-N7)):pi/2)).',1,size(MIn,2));


MOut=real(ifft(TF./abs(TF).*MOut));
