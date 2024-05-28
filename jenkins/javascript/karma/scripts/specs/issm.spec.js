//var json:any;
it("should load a fixture", function () {
    jasmine.getFixtures().fixturesPath = "base/Archives/"
    var f = readFixtures("Archive101.json");
    var json = JSON.parse(f);
    expect(json).toBeDefined();
    //console.log(json);
});

describe("test101", function() {
    it("contains test101", function() {
        console.log('creating model'); 
        var md = new model();

        console.log('meshing');
        triangle(md,square[0],40000); 

        console.log('parameterization');
        setmask(md,'all','');
        parameterize(md);
        setflowequation(md,'SSA','all');
        md.verbose.solution=2;  md.verbose.convergence=0;

        console.log("MESH");
        console.log(md.mesh.domaintype());

        md=solve(md,StressbalanceSolutionEnum(),'checkconsistency','no');

        console.log(md.results['StressbalanceSolution'][0]['Vel']);
    });
})

//describe("test105", function() {
    //it("contains test105", function() {
        ////Test Name: SquareShelfConstrainedMasstransp2d
        //var md = new model();
        //triangle(md,square[0],150000.);
        //setmask(md,'all','');
        //parameterize(md);
        //setflowequation(md,'SSA','all');
        ////md.cluster=generic('name',oshostname(),'np',3);
        //md=solve(md,MasstransportSolutionEnum());

        ////Fields and tolerances to track changes
        //field_names     =['Thickness'];
        //field_tolerances=[1e-13];
        //field_values=[
            //(md.results.MasstransportSolution[0].Thickness),
            //];
    //});
//})

//describe("test110", function() {
    //it("contains test110", function() {
        ////Test Name: SquareShelfConstrainedTranSSA2d
        //var md = new model();
        //triangle(md,square[0],150000.);
        //setmask(md,'all','');
        //parameterize(md);
        //setflowequation(md,'SSA','all');
        ////md.cluster=generic('name',oshostname(),'np',3);
        //md.trans.requested_outputs=['IceVolume'];

        //md=solve(md,TransientSolutionEnum());

        ////Fields and tolerances to track changes
        //field_names     =['Vx1','Vy1','Vel1','Pressure1','Bed1','Surface1','Thickness1','Volume1','Vx2','Vy2','Vel2','Pressure2','Bed2','Surface2','Thickness2','Volume2','Vx3','Vy3','Vel3','Pressure3','Bed3','Surface3','Thickness3','Volume3'];
        //field_tolerances=[1e-13,1e-13,1e-13,1e-13,1e-13,1e-13,1e-13,1e-13,
            //1e-13,1e-13,1e-13,1e-13,1e-13,1e-13,1e-13,1e-13,
        //1e-13,1e-13,1e-13,1e-13,1e-13,1e-13,1e-13,1e-13];
        //field_values=[
            //(md.results.TransientSolution[0](1).Vx),
            //(md.results.TransientSolution[0](1).Vy),
            //(md.results.TransientSolution[0](1).Vel),
            //(md.results.TransientSolution[0](1).Pressure),
            //(md.results.TransientSolution[0](1).Base),
            //(md.results.TransientSolution[0](1).Surface),
            //(md.results.TransientSolution[0](1).Thickness),
            //(md.results.TransientSolution[0](1).IceVolume),
            //(md.results.TransientSolution[0](2).Vx),
            //(md.results.TransientSolution[0](2).Vy),
            //(md.results.TransientSolution[0](2).Vel),
            //(md.results.TransientSolution[0](2).Pressure),
            //(md.results.TransientSolution[0](2).Base),
            //(md.results.TransientSolution[0](2).Surface),
            //(md.results.TransientSolution[0](2).Thickness),
            //(md.results.TransientSolution[0](2).IceVolume),
            //(md.results.TransientSolution[0](3).Vx),
            //(md.results.TransientSolution[0](3).Vy),
            //(md.results.TransientSolution[0](3).Vel),
            //(md.results.TransientSolution[0](3).Pressure),
            //(md.results.TransientSolution[0](3).Base),
            //(md.results.TransientSolution[0](3).Surface),
            //(md.results.TransientSolution[0](3).Thickness),
            //(md.results.TransientSolution[0](3).IceVolume),
        //];
    //});
//})

//describe("test112", function() {
    //it("contains test112", function() {
        ////Test Name: SquareShelfConstrainedSurfSlop2d
        //var md = new model();
        //triangle(md,square[0],150000.);
        //setmask(md,'all','');
        //parameterize(md);
        //setflowequation(md,'SSA','all');
        ////md.cluster=generic('name',oshostname(),'np',3);
        //md=solve(md,SurfaceSlopeSolutionEnum());

        ////Fields and tolerances to track changes
        //field_names     =['SurfaceSlopeX','SurfaceSlopeY'];
        //field_tolerances=[1e-13,1e-13];
        //field_values=[
            //(md.results.SurfaceSlopeSolution[0].SurfaceSlopeX),
            //(md.results.SurfaceSlopeSolution[0].SurfaceSlopeY),
        //];
    //});
//})

//describe("test114", function() {
    //it("contains test114", function() {
        ////Test Name: SquareShelfConstrainedBedSlop2d
        //var md = new model();
        //triangle(md,square[0],150000.);
        //setmask(md,'all','');
        //parameterize(md);
        //setflowequation(md,'SSA','all');
        ////md.cluster=generic('name',oshostname(),'np',3);
        //md=solve(md,BedSlopeSolutionEnum());

        ////Fields and tolerances to track changes
        //field_names     =['BedSlopeX','BedSlopeY'];
        //field_tolerances=[1e-13,1e-13];
        //field_values=[
            //(md.results.BedSlopeSolution[0].BedSlopeX),
            //(md.results.BedSlopeSolution[0].BedSlopeY),
        //];
    //});
//})

//describe("test201", function() {
    //it("contains test201", function() {
        ////Test Name: SquareShelfStressSSA2d
        //var md = new model();
        //triangle(md,square[0],150000.);
        //setmask(md,'all','');
        //parameterize(md);
        //setflowequation(md,'SSA','all');
        ////md.cluster=generic('name',oshostname(),'np',3);
        //md=solve(md,StressbalanceSolutionEnum());

        ////Fields and tolerances to track changes
        //field_names     =['Vx','Vy','Vel','Pressure'];
        //field_tolerances=[1e-13,1e-13,1e-13,1e-13];
        //field_values=[
            //(md.results.StressbalanceSolution[0].Vx),
            //(md.results.StressbalanceSolution[0].Vy),
            //(md.results.StressbalanceSolution[0].Vel),
            //(md.results.StressbalanceSolution[0].Pressure),
        //];
    //});
//})

//describe("test208", function() {
    //it("contains test208", function() {
        ////Test Name: SquareShelfTranSSA2d
        //var md = new model();
        //triangle(md,square[0],150000.);
        //setmask(md,'all','');
        //parameterize(md);
        //setflowequation(md,'SSA','all');
        ////md.cluster=generic('name',oshostname(),'np',3);
        //md.trans.requested_outputs=['default','FloatingArea','GroundedArea','TotalGroundedBmb','TotalFloatingBmb'];
        //for (var i = 0; i < md.basalforcings.floatingice_melting_rate.length; ++i) {
            //md.basalforcings.floatingice_melting_rate[i] = 1;
        //}
        //md=solve(md,TransientSolutionEnum());

        ////Fields and tolerances to track changes
        //field_names     =['Vx1','Vy1','Vel1','Pressure1','Bed1','Surface1','Thickness1','TotalGroundedBmb1','TotalFloatingBmb1',
            //'Vx2','Vy2','Vel2','Pressure2','Bed2','Surface2','Thickness2','TotalGroundedBmb2','TotalFloatingBmb2',
        //'Vx3','Vy3','Vel3','Pressure3','Bed3','Surface3','Thickness3','TotalGroundedBmb3','TotalFloatingBmb3'];
        //field_tolerances=[1e-13,1e-13,1e-13,1e-13,1e-13,1e-13,1e-13,1e-13,1e-13,
            //1e-13,1e-13,1e-13,1e-13,1e-13,1e-13,1e-13,1e-13,1e-13,
        //1e-13,1e-13,1e-13,1e-13,1e-13,1e-13,1e-13,1e-13,1e-13];
        //field_values=[
            //(md.results.TransientSolution[0](1).Vx),
            //(md.results.TransientSolution[0](1).Vy),
            //(md.results.TransientSolution[0](1).Vel),
            //(md.results.TransientSolution[0](1).Pressure),
            //(md.results.TransientSolution[0](1).Base),
            //(md.results.TransientSolution[0](1).Surface),
            //(md.results.TransientSolution[0](1).Thickness),
            //(md.results.TransientSolution[0](1).TotalGroundedBmb),
            //(md.results.TransientSolution[0](1).TotalFloatingBmb),
            //(md.results.TransientSolution[0](2).Vx),
            //(md.results.TransientSolution[0](2).Vy),
            //(md.results.TransientSolution[0](2).Vel),
            //(md.results.TransientSolution[0](2).Pressure),
            //(md.results.TransientSolution[0](2).Base),
            //(md.results.TransientSolution[0](2).Surface),
            //(md.results.TransientSolution[0](2).Thickness),
            //(md.results.TransientSolution[0](2).TotalGroundedBmb),
            //(md.results.TransientSolution[0](2).TotalFloatingBmb),
            //(md.results.TransientSolution[0](3).Vx),
            //(md.results.TransientSolution[0](3).Vy),
            //(md.results.TransientSolution[0](3).Vel),
            //(md.results.TransientSolution[0](3).Pressure),
            //(md.results.TransientSolution[0](3).Base),
            //(md.results.TransientSolution[0](3).Surface),
            //(md.results.TransientSolution[0](3).Thickness),
            //(md.results.TransientSolution[0](3).TotalGroundedBmb),
            //(md.results.TransientSolution[0](3).TotalFloatingBmb),
        //];
    //});
//})
