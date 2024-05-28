%Quick documentation for ISSM

%First get ISSM tier: 
ISSM_DIR=issmdir();

disp('  A comprehensive documentation is available on http://issm.jpl.nasa.gov');
disp('  Example: how to create a square ice shelf');
disp(['       go to ',ISSM_DIR,'/examples/SquareIceShelf']);
disp(sprintf('%-63s %s','       md=model;','%creates a new empty model structure'));
disp(sprintf('%-63s %s','       md=triangle(md,''DomainOutline.exp'',50000);','%creates a mesh of the domain outline with a resolution of 50000 m'));
disp(sprintf('%-63s %s','       md=setmask(md,''all'','''');','%defines the glacier system as an ice shelf (no island)'));
disp(sprintf('%-63s %s','       md=parameterize(md,''Square.par'');','%fills all the other fields of the model'));
disp(sprintf('%-63s %s','       md=setflowequation(md,''SSA'',''all'');','%defines all elements as SSA''s SSA'));
disp(sprintf('%-63s %s','       md=solve(md,''Stressbalance'');','%solve for stress balance'));
disp(sprintf('%-63s %s','       plotmodel(md,''data'',md.results.StressbalanceSolution.Vel);','%displays the velocity (type plotdoc for plotmodel help)'));
