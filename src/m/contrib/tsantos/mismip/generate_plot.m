%L0 and L1 {{{
md0			= loadmodel('./L0_viscous/Transient_Steadystate_Level_0.mat');
md1			= loadmodel('./L1/L1_viscous/Transient_steadystate.mat');
md_L1_R5		= loadmodel('./R5/L1_viscous/Transient_steadystate.mat');
md_L1_R15	= loadmodel('./R15/L1_viscous/Transient_steadystate.mat');
md_L1_R30	= loadmodel('./R30/L1_viscous/Transient_steadystate.mat');
md_L1_R40	= loadmodel('./R40/L1_viscous/Transient_steadystate.mat');

[ga0 iv0 ivaf0 GLy400 ne0 t0] = ice_evolution(md0);
[gaL1 ivL1 ivafL1 GLy40L1 neL1 tL1] = ice_evolution(md1);
[ga1 iv1 ivaf1 GLy401 ne1 t1] = ice_evolution(md_L1_R5);
[ga2 iv2 ivaf2 GLy402 ne2 t2] = ice_evolution(md_L1_R15);
[ga3 iv3 ivaf3 GLy403 ne3 t3] = ice_evolution(md_L1_R30);
[ga4 iv4 ivaf4 GLy404 ne4 t4] = ice_evolution(md_L1_R40);
% }}}
%L2 {{{
md_L2_R5		= loadmodel('./R5/L2_viscous/Transient_steadystate.mat');
md_L2_R15	= loadmodel('./R15/L2_viscous/Transient_steadystate.mat');
md_L2_R20	= loadmodel('./R20/L2_viscous/Transient_steadystate.mat');
md_L2_R30	= loadmodel('./R30/L2_viscous/Transient_steadystate.mat');

[ga5 iv5 ivaf5 GLy405 ne5 t5] = ice_evolution(md_L2_R5);
[ga6 iv6 ivaf6 GLy406 ne6 t6] = ice_evolution(md_L2_R15);
[ga7 iv7 ivaf7 GLy407 ne7 t7] = ice_evolution(md_L2_R20);
[ga8 iv8 ivaf8 GLy408 ne8 t8] = ice_evolution(md_L2_R30);
% }}}
%L3 {{{
md_L3_R30	= loadmodel('./L3_viscous/Transient_steadystate.mat');
[ga9 iv9 ivaf9 GLy409 ne9 t9] = ice_evolution(md_L3_R30);
% }}}

%L4 {{{
md_L4_R30	= loadmodel('./L4_viscous/Transient_steadystate.mat');
[ga10 iv10 ivaf10 GLy4010 ne10 t10] = ice_evolution(md_L4_R30);
% }}}

% Scaling {{{
t0=t0/1000;
tL1=tL1/1000;
t1=t1/1000;
t2=t2/1000;
t3=t3/1000;
t4=t4/1000;
t5=t5/1000;
t6=t6/1000;
t7=t7/1000;
t8=t8/1000;
t9=t9/1000;
t10=t10/1000;

GLy400=GLy400/1000;
GLy40L1=GLy40L1/1000;
GLy401=GLy401/1000;
GLy402=GLy402/1000;
GLy403=GLy403/1000;
GLy404=GLy404/1000;
GLy405=GLy405/1000;
GLy406=GLy406/1000;
GLy407=GLy407/1000;
GLy408=GLy408/1000;
GLy409=GLy409/1000;
GLy4010=GLy4010/1000;

ivaf0=ivaf0/10^12;
ivafL1=ivafL1/10^12;
ivaf1=ivaf1/10^12;
ivaf2=ivaf2/10^12;
ivaf3=ivaf3/10^12;
ivaf4=ivaf4/10^12;
ivaf5=ivaf5/10^12;
ivaf6=ivaf6/10^12;
ivaf7=ivaf7/10^12;
ivaf8=ivaf8/10^12;
ivaf9=ivaf9/10^12;
ivaf10=ivaf10/10^12;
% }}}

figure(1) % {{{
hold on
plot(t0,GLy400,'*',tL1,GLy40L1,'r+',t1,GLy401,'b',t2,GLy402,'g',t3,GLy403,'m',t4,GLy404,'r','LineWidth',1.5)
axis([0 25 365 455])
legend('L0','L1','L1-R5','L1-R15','L1-R30','L1-R40','Location','southeast')
title('\fontsize{16} Mismip+ steady state - Rmax analysis')
ylabel('GL @ y=40km (km)')
xlabel('t (kyr)')
hold off
% }}}

figure(2) % {{{
hold on
plot(t0,GLy400,'*',t5,GLy405,'b',t6,GLy406,'g',t7,GLy407,'m',t8,GLy408,'r','LineWidth',1.5)
axis([0 25 365 455])
legend('L0','L2-R5','L2-R15','L2-R20','L2-R30','Location','southeast')
title('\fontsize{16} Mismip+ steady state - Rmax analysis')
ylabel('GL @ y=40km (km)')
xlabel('t (kyr)')
hold off
% }}}

figure(3) % {{{
res	= [2 1];
res30 = [2 1 0.5 0.25];
r5		= [GLy401(end) GLy405(end)];
r15	= [GLy402(end) GLy406(end)];
r30	= [GLy403(end) GLy408(end) GLy409(end) GLy4010(end)];
hold on
plot(log(4),GLy400(end),'^','MarkerSize',7,'MarkerFaceColor','b','MarkerEdgeColor','k');
plot(log(2),GLy40L1(end),'^','MarkerSize',7,'MarkerFaceColor','k','MarkerEdgeColor','k');
plot(log(res),r5,'^','MarkerSize',7,'MarkerFaceColor','g','MarkerEdgeColor','k');
plot(log(res),r15,'^','MarkerSize',7,'MarkerFaceColor','y','MarkerEdgeColor','k');
plot(log(res30),r30,'^','MarkerSize',7,'MarkerFaceColor','r','MarkerEdgeColor','k');
plot(log(res(1)),GLy404(end),'^','MarkerSize',7,'MarkerFaceColor','m','MarkerEdgeColor','k');
plot(log(res(2)),GLy407(end),'^','MarkerSize',7,'MarkerFaceColor','c','MarkerEdgeColor','k');
axis([log(0.2) log(5) 437 458])
legend('L0','L1','R5','R15','R30','L1-R40','L2-R20','Location','northeast')
title('\fontsize{15} Mismip+ steady state: GLpos X Rmax')
ylabel('GL @ y=40km (km)')
xlabel('Resolution (km)')
xtickslabel = [0.25 0.5 1 2 4];
xticks = log(xtickslabel);
ytickslabel = [440 445 450 455];
yticks = ytickslabel;
ax = gca;
set(ax,'XTick',xticks);
set(ax,'XTickLabel',xtickslabel);
set(ax,'YTick',yticks);
set(ax,'YTickLabel',ytickslabel);
box on
%set(gca,'xscale','log')
hold off
% }}}

figure(4) % {{{
res	= [2 1];
res30	= [2 1 0.5 0.25];
r5		= [ivaf1(end) ivaf5(end)];
r15	= [ivaf2(end) ivaf6(end)];
r30	= [ivaf3(end) ivaf8(end) ivaf9(end) ivaf10(end)];
hold on
plot(log(4),ivaf0(end),'^','MarkerSize',7,'MarkerFaceColor','b','MarkerEdgeColor','k');
plot(log(2),ivafL1(end),'^','MarkerSize',7,'MarkerFaceColor','k','MarkerEdgeColor','k');
plot(log(res),r5,'^','MarkerSize',7,'MarkerFaceColor','g','MarkerEdgeColor','k');
plot(log(res),r15,'^','MarkerSize',7,'MarkerFaceColor','y','MarkerEdgeColor','k');
plot(log(res30),r30,'^','MarkerSize',7,'MarkerFaceColor','r','MarkerEdgeColor','k');
plot(log(res(1)),ivaf4(end),'^','MarkerSize',7,'MarkerFaceColor','m','MarkerEdgeColor','k');
plot(log(res(2)),ivaf7(end),'^','MarkerSize',7,'MarkerFaceColor','c','MarkerEdgeColor','k');
axis([log(0.2) log(5) 34.3 37.2])
legend('L0','L1','R5','R15','R30','L1-R40','L2-R20','Location','northeast')
title('\fontsize{15} Mismip+ steady state: IVAF X Rmax')
ylabel('Ice volume above floatation (1000 x km3)')
xlabel('Resolution (km)')
xtickslabel = [0.25 0.5 1 2 4];
xticks = log(xtickslabel);
ytickslabel = [34.5 35.0 35.5 36.0 36.5 37.0];
yticks = ytickslabel;
ax = gca;
set(ax,'XTick',xticks);
set(ax,'XTickLabel',xtickslabel);
set(ax,'YTick',yticks);
set(ax,'YTickLabel',ytickslabel);
box on
%set(gca,'xscale','log')
hold off
% }}}

figure(5) % {{{
hold on
plot(t0,GLy400,'*',t3,GLy403,'b',t8,GLy408,'m',t9,GLy409,'r',t10,GLy4010,'g','LineWidth',1.5)
axis([0 25 365 460])
legend('L0','L1-R30','L2-R30','L3-R30','L4-R30','Location','southeast')
title('\fontsize{16} Mismip+ steady state - Rmax analysis')
ylabel('GL @ y=40km (km)')
xlabel('t (kyr)')
hold off
% }}}
