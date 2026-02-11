% read test4003 MItgcm output
ISSM_DIR=getenv('ISSM_DIR');
eval (['cd ' ISSM_DIR '/test/NightlyRun/run'])
fnm=dir('surfDiag.*.data');
fld={'SHIfwFlx (kg/m^2/s) Ice shelf fresh water flux (positive upward)', ...
     'SHIhtFlx (W/m^2   ) Ice shelf heat flux  (positive upward)', ...
     'SHIgammT (m/s     ) Ice shelf exchange coefficient for theta', ...
     'SHIgammS (m/s     ) Ice shelf exchange coefficient for salt', ...
     'SHI_mass (kg/m^2  ) dynamic ice shelf mass for surface load anomaly', ...
     'SHIRshel (m       ) depth of shelfice', ...
     'SI_Uvel  (m/a     ) Ice stream x-velocity', ...
     'SI_Vvel  (m/a     ) Ice stream y-velocity', ...
     'SI_Thick (m       ) Ice stream thickness', ...
     'SI_hmask (none    ) Ice stream thickness mask', ...
     'SI_float (none    ) Ice stream grounding ind', ...
     'SHIuStar (m/s     ) Friction velocity at bottom of ice shelf'}
nx=3; ny=200; nf=length(fld); nt=length(fnm);
Diag=zeros(nx,ny,nf,nt);
for t=1:nt
    Diag(:,:,:,t)=readbin(fnm(t).name,[nx ny nf]);
end

% plot output
mkdir figs
orient tall
wysiwyg
colormap(cmap)
for t=1:nt
    clf
    for f=1:nf
        subplot(nf,1,f)
        if t==1
            mypcolor(Diag(:,:,f,t));
        else
            mypcolor(Diag(:,:,f,t)-Diag(:,:,f,1));
        end
        colorbar
        title(fld{f},'Interpreter','none')
        set(gca,'YTickLabel',[])
        if f<nf
            set(gca,'XTickLabel',[])
        end
    end
    eval(['print -djpeg figs/surfDiag' myint2str(t)])
end
