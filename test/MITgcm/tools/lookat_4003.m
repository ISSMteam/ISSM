pn='~/mitgcm/slr/components/issm/trunk-jpl/test/NightlyRun/';
p1=[pn 'RunUncoupled/'];
p2=[pn 'run/'];

for ts=30:34
    v1=readbin([p2 'R_shelfIce1_' myint2str(ts,10) '.data'],[3 200]);
    v2=readbin([p2 'R_shelfIce2_' myint2str(ts,10) '.data'],[3 200]);
    u=readbin([p2 'U.' myint2str(ts,10) '.data'],[3 200 90]);
    v=readbin([p2 'V.' myint2str(ts,10) '.data'],[3 200 90]);
    w=readbin([p2 'W.' myint2str(ts,10) '.data'],[3 200 90]);
    clf, subplot(511), plot(v1(2,:)), title(ts)
    subplot(512), plot(v2(2,:)-v1(2,:)), title('draft change')
    subplot(513), mypcolor(1:200,-1:-1:-90,squeeze(u(2,:,:))'); title('U'), colorbar
    subplot(514), mypcolor(1:200,-1:-1:-90,squeeze(v(2,:,:))'); title('V'), colorbar
    subplot(515), mypcolor(1:200,-1:-1:-90,squeeze(w(2,:,:))'); title('W'), colorbar
    pause
end


    clf
    subplot(311), mypcolor(v1), title(ts)
    subplot(312), mypcolor(v2), title(ts)
    subplot(313), mypcolor(v2-v1)
    
for ts=1:10  
    v1=readbin([p2 'SHICE_fwFlux.' myint2str(ts,10) '.data'],[3 200]);
    clf
    subplot(311), mypcolor(v1), title(ts), colorbar
pause
end


    v2=readbin([p2 'SHICE_fwFlux.' myint2str(ts,10) '.data'],[3 200]);


clf
for ts=5, disp(ts)
    T=readbin([p2 'T.' myint2str(ts,10) '.data'],[3 200 90]);
    S=readbin([p2 'S.' myint2str(ts,10) '.data'],[3 200 90]);
    U=readbin([p2 'U.' myint2str(ts,10) '.data'],[3 200 90]);
    V=readbin([p2 'V.' myint2str(ts,10) '.data'],[3 200 90]);
    W=readbin([p2 'W.' myint2str(ts,10) '.data'],[3 200 90]);
    for k=1:90, disp(k)
        clf
        subplot(321), plot(S(2,:,k)), title('S')
        subplot(322), plot(T(2,:,k)), title('T')
        subplot(323), plot(U(2,:,k)), title('U')
        subplot(324), plot(V(2,:,k)), title('V')
        subplot(325), plot(W(2,:,k)), title('W')
        pause
    end
end


clf
for ts=1:8, disp(ts)
    T=readbin([p2 'T.' myint2str(ts,10) '.data'],[3 200 90]);
    S=readbin([p2 'S.' myint2str(ts,10) '.data'],[3 200 90]);
    U=readbin([p2 'U.' myint2str(ts,10) '.data'],[3 200 90]);
    V=readbin([p2 'V.' myint2str(ts,10) '.data'],[3 200 90]);
    W=readbin([p2 'W.' myint2str(ts,10) '.data'],[3 200 90]);
    clf
    subplot(321), mypcolor(squeeze(S(2,:,:))'), title('S'), colorbar
    subplot(322), mypcolor(squeeze(T(2,:,:))'), title('T'), colorbar
    subplot(323), mypcolor(squeeze(U(2,:,:))'), title('U'), colorbar
    subplot(324), mypcolor(squeeze(V(2,:,:))'), title('V'), colorbar
    subplot(325), mypcolor(squeeze(W(2,:,:))'), title('W'), colorbar
    pause
end






ts=8;
v1=readbin([p2 'R_shelfIce1_' myint2str(ts,10) '.data'],[3 200]);
v2=readbin([p2 'R_shelfIce2_' myint2str(ts,10) '.data'],[3 200]);
clf
subplot(311), mypcolor(v1); title('R_shelfIce1'), colorbar
subplot(312), mypcolor(v2); title('R_shelfIce2'), colorbar
subplot(313), mypcolor(v2-v1); title('diff'), colorbar

figure(2)
clf
plot(1:200,v1(2,:),'o-',1:200,v2(2,:),'o-',1:200,v2(2,:)-v1(2,:),'o-')

ts=0;
fld='R_shelfIce1_';
v1=readbin([p2 fld myint2str(ts,10) '.data'],[3 200]);
fld='R_shelfIce2_';
for ts=0:8:184
    v2=readbin([p2 fld myint2str(ts,10) '.data'],[3 200]);
    clf
    subplot(311), mypcolor(v1); title(ts-8), colorbar
    subplot(312), mypcolor(v2); title(ts), colorbar
    subplot(313), mypcolor(v2-v1); title('diff'), colorbar
    pause
    v1=v2;
end

fld='surfDiag';
ts=2;
v1=rdmds([p1 fld],ts);
v2=rdmds([p2 fld],ts);
clf
fld={'ETAN','RSURF','oceTAUX','oceTAUY','oceQnet','oceFWflx', ...
     'MXLDEPTH','SHIfwFlx','SHIhtFlx','SHIgammT','SHIgammS', ...
     'SHI_mass','SHIuStar'};
for i=1:length(fld)
    subplot(311), mypcolor(v1(:,:,i)); title(fld{i}), colorbar
    subplot(312), mypcolor(v2(:,:,i)); title('coupled'), colorbar
    subplot(313), mypcolor(v2(:,:,i)-v1(:,:,i)); title('coupled-uncoupled'), colorbar
    pause
end

fld='Eta';
ts=1;
v1=rdmds([p1 fld],ts);
v2=rdmds([p2 fld],ts);
clf
subplot(311), mypcolor(v1); title(fld), colorbar
subplot(312), mypcolor(v2); title('coupled'), colorbar
subplot(313), mypcolor(v2-v1); title('coupled-uncoupled'), colorbar
