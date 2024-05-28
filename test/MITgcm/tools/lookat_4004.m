clear, close all
pn='~/mitgcm/slr/components/issm/trunk-jpl/test/NightlyRun/run/';
figure(1), clf, orient tall, wysiwyg
for ts=0:4
    r1=readbin([pn 'R_shelfIce1_' myint2str(ts+1,10) '.data'],[3 200]);
    r2=readbin([pn 'R_shelfIce2_' myint2str(ts,10) '.data'],[3 200]);
    subplot(5,1,ts+1), mypcolor(r2-r1); colorbar
end

for ts=1:5
    figure(ts+1), clf, orient tall, wysiwyg
    e=readbin([pn 'Eta.' myint2str(ts,10) '.data'],[3 200]);
    p=readbin([pn 'PHL.' myint2str(ts,10) '.data'],[3 200]);
    f=readbin([pn 'SHICE_fwFlux.' myint2str(ts,10) '.data'],[3 200]);
    h=readbin([pn 'SHICE_heatFlux.' myint2str(ts,10) '.data'],[3 200]);
    ph=readbin([pn 'PH.' myint2str(ts,10) '.data'],[3 200 90]);
    clf, subplot(511), mypcolor(e); title(['Eta @ ts=' int2str(ts)]), colorbar
    subplot(512), mypcolor(p); title('PHL'), colorbar
    subplot(513), mypcolor(f); title('SHICE fwFlux'), colorbar
    subplot(514), mypcolor(h); title('SHICE heatFlux'), colorbar
    subplot(515), mypcolor(1:200,-1:-1:-90,squeeze(ph(2,:,:))'); title('PH'), colorbar
end

for ts=0:5
    figure(ts+7), clf, orient tall, wysiwyg
    r2=readbin([pn 'R_shelfIce2_' myint2str(ts,10) '.data'],[3 200]);
    if ts<5
        r1=readbin([pn 'R_shelfIce1_' myint2str(ts+1,10) '.data'],[3 200]);
    end
    s =readbin([pn 'S.' myint2str(ts,10) '.data'],[3 200 90]);
    t =readbin([pn 'T.' myint2str(ts,10) '.data'],[3 200 90]);
    u =readbin([pn 'U.' myint2str(ts,10) '.data'],[3 200 90]);
    v =readbin([pn 'V.' myint2str(ts,10) '.data'],[3 200 90]);
    w =readbin([pn 'W.' myint2str(ts,10) '.data'],[3 200 90]);
    in=find(~s); s(in)=nan; t(in)=nan;
    clf,
    subplot(711), mypcolor(r2); title(['draft @ ts=' int2str(ts)]), colorbar
    subplot(712), mypcolor(r2-r1); title('draft change'), colorbar
    subplot(713), pcolorcen(1:200,-1:-1:-90,squeeze(s(2,:,:))'); title('S'), colorbar
    subplot(714), pcolorcen(1:200,-1:-1:-90,squeeze(t(2,:,:))'); title('T'), colorbar
    subplot(715), pcolorcen(1:200,-1:-1:-90,squeeze(u(2,:,:))'); title('U'), colorbar
    subplot(716), pcolorcen(1:200,-1:-1:-90,squeeze(v(2,:,:))'); title('V'), colorbar
    subplot(717), pcolorcen(1:200,-1:-1:-90,squeeze(w(2,:,:))'); title('W'), colorbar
end
