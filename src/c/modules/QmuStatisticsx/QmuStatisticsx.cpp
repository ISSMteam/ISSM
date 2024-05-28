/*!\file:  QmuStatisticsx routines
 */ 
/*includes and prototypes:*/
#include <sys/stat.h>
#include "./QmuStatisticsx.h"
#include "../OutputResultsx/OutputResultsx.h"

int readdata(IssmDouble** pdoublemat, int* pdoublematsize, IssmDouble* pdouble, FILE* fid,char* field,int step){ /*{{{*/

	int length;
	char fieldname[1000];
	int   fieldname_size;
	IssmDouble   rtime;
	int          rstep;
	int M,N;

	//fields that we retrive: 
	IssmDouble  dfield; 
	char*       sfield    = NULL;
	IssmDouble* dmatfield = NULL; 
	int*        imatfield = NULL; 

	//type of the returned field: 
	int type;
	int found=0;

	while(1){

		size_t ret_code = fread(&fieldname_size, sizeof(int), 1, fid); 
		if(ret_code != 1) break; //we are done.

		fread(fieldname, sizeof(char), fieldname_size, fid); 
		//_printf0_("fieldname: " << fieldname << "\n");

		fread(&rtime, sizeof(IssmDouble), 1, fid); 
		fread(&rstep, sizeof(int), 1, fid); 

		//check on field: 
		if ((step==rstep) && (strcmp(field,fieldname)==0)){

			//ok, go read the result really: 
			fread(&type,sizeof(int),1,fid);
			fread(&M,sizeof(int),1,fid);
			if (type==1){
				fread(&dfield,sizeof(IssmDouble),1,fid);
			}
			else if (type==2){
				fread(&M,sizeof(int),1,fid);
				sfield=xNew<char>(M);
				fread(sfield,sizeof(char),M,fid);
			}
			else if (type==3){
				fread(&N,sizeof(int),1,fid);
				dmatfield=xNew<IssmDouble>(M*N);
				fread(dmatfield,sizeof(IssmDouble),M*N,fid);
			}
			else if (type==4){
				fread(&N,sizeof(int),1,fid);
				imatfield=xNew<int>(M*N);
				fread(imatfield,sizeof(int),M*N,fid);
			}
			else _error_("cannot read data of type " << type << "\n");
			found=1;
			break;
		}
		else{
			//just skim to next results.
			fread(&type,sizeof(int),1,fid);
			fread(&M,sizeof(int),1,fid);
			if (type==1){
				fseek(fid,sizeof(IssmDouble),SEEK_CUR);
			}
			else if(type==2){
				fseek(fid,M*sizeof(char),SEEK_CUR);
			}
			else if(type==3){
				fread(&N,sizeof(int),1,fid);
				fseek(fid,M*N*sizeof(IssmDouble),SEEK_CUR);
			}
			else if(type==4){
				fread(&N,sizeof(int),1,fid);
				fseek(fid,M*N*sizeof(int),SEEK_CUR);
			}
			else _error_("cannot read data of type " << type << "\n");
		}
	}
	if(found==0)_error_("could not find " << field << " at step " << step  << "\n");

	/*assign output pointers:*/
	*pdoublemat=dmatfield;
	*pdoublematsize=M*N;
	*pdouble=dfield;

	/*return:*/
	return type;

}
/*}}}*/
int ComputeHistogram(Parameters* parameters,Results* results,int color, ISSM_MPI_Comm statcomm){  /*{{{*/

	int nsamples; 
	char* directory=NULL;
	char* model=NULL;
	char** fields=NULL;
	int* steps=NULL;
	int nsteps;
	int nfields;
	int nbins;
	int range,lower_row,upper_row;
	int nfilesperdirectory;

	/*intermediary:*/
	IssmDouble* doublemat=NULL;
	int         doublematsize;
	IssmDouble scalar;

	/*computation of average and variance itself:*/
	IssmDouble** maxxs = NULL;
	IssmDouble** minxs = NULL;
	int*         xtype=NULL;
	int*         xsize=NULL;

	IssmDouble** maxmeans=NULL;
	IssmDouble** minmeans=NULL;
	int*         meanxtype=NULL;
	int*         meanxsize=NULL;

	/*only work on the statistical communicator: */
	if (color==MPI_UNDEFINED)return 0;

	/*Retrieve parameters:*/
	parameters->FindParam(&nfilesperdirectory,QmuNfilesPerDirectoryEnum);
	parameters->FindParam(&nsamples,QmuNsampleEnum);
	parameters->FindParam(&directory,DirectoryNameEnum);
	parameters->FindParam(&model,InputFileNameEnum);
	parameters->FindParam(&fields,&nfields,FieldsEnum);
	parameters->FindParam(&steps,&nsteps,StepsEnum);
	parameters->FindParam(&nbins,NbinsEnum);

	/*Get rank from the stat comm communicator:*/
	IssmComm::SetComm(statcomm);
	int my_rank=IssmComm::GetRank();

	/*Open files and read them complelety, in a distributed way:*/
	range=DetermineLocalSize(nsamples,IssmComm::GetComm());
	GetOwnershipBoundariesFromRange(&lower_row,&upper_row,range,IssmComm::GetComm());

	/*Initialize arrays:*/
	maxmeans=xNew<IssmDouble*>(nfields);
	minmeans=xNew<IssmDouble*>(nfields);
	meanxtype=xNew<int>(nfields);
	meanxsize=xNew<int>(nfields);

	maxxs=xNew<IssmDouble*>(nfields*nsteps);
	minxs=xNew<IssmDouble*>(nfields*nsteps);
	xtype=xNew<int>(nfields*nsteps);
	xsize=xNew<int>(nfields*nsteps);

	/*Start opening files:*/
	for(int i=(lower_row+1);i<=upper_row;i++){
		_printf0_("reading file #: " << i << "\n");
		char file[1000];
		long int  length;
		char* buffer=NULL;

		/*string:*/
		sprintf(file,"%s/%i/%s.outbin.%i",directory,my_rank+1,model,i);

		/*open file: */
		_printf0_("    opening file: " << file << "\n");
		FILE* fid=fopen(file,"rb");
		if(fid==NULL)_error_("could not open file: " << file << "\n");

		/*figure out size of file, and read the whole thing:*/
		_printf0_("    reading file:\n");
		fseek(fid, 0, SEEK_END);
		length = ftell (fid);
		fseek(fid, 0, SEEK_SET);
		buffer = xNew<char>(length);
		fread(buffer, sizeof(char), length, fid);

		/*close file:*/
		fclose(fid);

		/*create a memory stream with this buffer:*/
		_printf0_("    processing file:\n");
		fid=fmemopen(buffer, length, "rb");

		/*start reading data from the buffer directly:*/
		for (int f=0;f<nfields;f++){
			char* field=fields[f];
			fseek(fid,0,SEEK_SET);
			for (int j=0;j<nsteps;j++){
				int counter=f*nsteps+j;
				xtype[counter]=readdata(&doublemat, &doublematsize, &scalar, fid,field,steps[j]);
				if(i==(lower_row+1)){
					if(xtype[counter]==1){
						maxxs[counter]=xNew<IssmDouble>(1); 
						minxs[counter]=xNew<IssmDouble>(1); 
						*maxxs[counter]=scalar;
						*minxs[counter]=scalar;
						xsize[counter]=1;
					}
					else if (xtype[counter]==3){
						maxxs[counter]=xNew<IssmDouble>(doublematsize); 
						xMemCpy<IssmDouble>(maxxs[counter],doublemat,doublematsize);
						minxs[counter]=xNew<IssmDouble>(doublematsize); 
						xMemCpy<IssmDouble>(minxs[counter],doublemat,doublematsize);
						xsize[counter]=doublematsize;
						xDelete<IssmDouble>(doublemat);
					}
					else _error_("cannot carry out statistics on type " << xtype[counter]); 
				}
				else{
					if(xtype[counter]==1){
						*maxxs[counter]=max(*maxxs[counter],scalar);
						*minxs[counter]=min(*minxs[counter],scalar);
					}
					else if (xtype[counter]==3){
						IssmDouble* newmax=maxxs[counter];
						IssmDouble* newmin=minxs[counter];
						for(int k=0;k<doublematsize;k++){
							if(doublemat[k]>newmax[k])newmax[k]=doublemat[k];
							if(doublemat[k]<newmin[k])newmin[k]=doublemat[k];
						}
						xDelete<IssmDouble>(doublemat);
					}
					else _error_("cannot carry out statistics on type " << xtype[counter]); 
				}
			}
		}
		_printf0_("    average in time:\n");

		/*Deal with average in time: */
		for (int f=0;f<nfields;f++){
			fseek(fid,0,SEEK_SET);
			char* field=fields[f];
			meanxtype[f]=readdata(&doublemat, &doublematsize, &scalar, fid,field,steps[0]);

			if(meanxtype[f]==1){
				meanxsize[f]=1;
				IssmDouble timemean=0;
				fseek(fid,0,SEEK_SET);
				for (int j=0;j<nsteps;j++){
					readdata(&doublemat, &doublematsize, &scalar, fid,field,steps[j]);
					timemean+=scalar/nsteps;
				}

				/*Figure out max and min of time means: */
				if(i==(lower_row+1)){
					maxmeans[f]=xNewZeroInit<IssmDouble>(1); 
					minmeans[f]=xNewZeroInit<IssmDouble>(1); 
					*maxmeans[f]=timemean;
					*minmeans[f]=timemean;
				}
				else{
					*maxmeans[f]=max(*maxmeans[f],timemean);
					*minmeans[f]=min(*minmeans[f],timemean);
				}
			}
			else{
				meanxsize[f]=doublematsize;
				fseek(fid,0,SEEK_SET);
				IssmDouble* timemean=xNewZeroInit<IssmDouble>(doublematsize);
				for (int j=0;j<nsteps;j++){
					readdata(&doublemat, &doublematsize, &scalar, fid,field,steps[j]);
					for (int k=0;k<doublematsize;k++){
						timemean[k]+=doublemat[k]/nsteps;
					}
					xDelete<IssmDouble>(doublemat);
				}

				if(i==(lower_row+1)){
					maxmeans[f]=xNew<IssmDouble>(doublematsize);
					xMemCpy<IssmDouble>(maxmeans[f],timemean,doublematsize);
					minmeans[f]=xNew<IssmDouble>(doublematsize);
					xMemCpy<IssmDouble>(minmeans[f],timemean,doublematsize);
				}
				else{
					IssmDouble* maxx=maxmeans[f];
					IssmDouble* minx=minmeans[f];

					for(int k=0;k<doublematsize;k++){
						maxx[k]=max(maxx[k],timemean[k]);
						minx[k]=min(minx[k],timemean[k]);
					}
					maxmeans[f]=maxx;
					minmeans[f]=minx;
				}
			}
		}
		fclose(fid);

		/*delete buffer:*/
		xDelete<char>(buffer);
	}
	ISSM_MPI_Barrier(IssmComm::GetComm());
	_printf0_("Done reading files, now computing min and max.\n"); 

	/*We have agregated minx and max across the cluster, now gather across the cluster onto
	 *cpu0 and then compute statistics:*/
	for (int f=0;f<nfields;f++){
		int counter0=f*nsteps+0;
		if (xtype[counter0]==1){ /*deal with scalars {{{*/
			for (int j=0;j<nsteps;j++){
				int counter=f*nsteps+j;

				/*we are broadcasting doubles:*/
				IssmDouble maxscalar=*maxxs[counter];
				IssmDouble minscalar=*minxs[counter];
				IssmDouble allmaxscalar;
				IssmDouble allminscalar;
				IssmDouble sumscalar_alltimes=0;

				ISSM_MPI_Allreduce(&maxscalar,&allmaxscalar,1,ISSM_MPI_PDOUBLE,ISSM_MPI_MAX,IssmComm::GetComm());
				ISSM_MPI_Allreduce(&minscalar,&allminscalar,1,ISSM_MPI_PDOUBLE,ISSM_MPI_MIN,IssmComm::GetComm());

				/*Store broadcasted value for later computation of histograms:*/
				*maxxs[counter]=allmaxscalar;
				*minxs[counter]=allminscalar;

			}
		} /*}}}*/
		else{ /*deal with arrays:{{{*/

			int size=xsize[counter0];
			for (int j=0;j<nsteps;j++){
				int counter=f*nsteps+j;

				/*we are broadcasting double arrays:*/
				IssmDouble* maxx=maxxs[counter];
				IssmDouble* minx=minxs[counter];

				IssmDouble*  allmax=xNew<IssmDouble>(size);
				IssmDouble*  allmin=xNew<IssmDouble>(size);

				ISSM_MPI_Allreduce(maxx,allmax,size,ISSM_MPI_PDOUBLE,ISSM_MPI_MAX,IssmComm::GetComm());
				ISSM_MPI_Allreduce(minx,allmin,size,ISSM_MPI_PDOUBLE,ISSM_MPI_MIN,IssmComm::GetComm());

				/*Store broadcasted value for later computation of histograms:*/
				maxxs[counter]=allmax;
				minxs[counter]=allmin;
			}
		} /*}}}*/
	}

	/*Now do the same for the time mean fields:*/
	for (int f=0;f<nfields;f++){
		if (meanxtype[f]==1){ /*deal with scalars {{{*/

			/*we are broadcasting doubles:*/
			IssmDouble maxscalar=*maxmeans[f];
			IssmDouble minscalar=*minmeans[f];
			IssmDouble allmaxscalar;
			IssmDouble allminscalar;

			ISSM_MPI_Allreduce(&maxscalar,&allmaxscalar,1,ISSM_MPI_PDOUBLE,ISSM_MPI_MAX,IssmComm::GetComm());
			ISSM_MPI_Allreduce(&minscalar,&allminscalar,1,ISSM_MPI_PDOUBLE,ISSM_MPI_MIN,IssmComm::GetComm());

			/*Store for later use in histogram computation:*/
			*maxmeans[f]=allmaxscalar;
			*minmeans[f]=allminscalar;

		} /*}}}*/
		else{ /*deal with arrays:{{{*/

			int size=meanxsize[f];

			/*we are broadcasting double arrays:*/
			IssmDouble* maxx=maxmeans[f];
			IssmDouble* minx=minmeans[f];

			IssmDouble*  allmax=xNew<IssmDouble>(size);
			IssmDouble*  allmin=xNew<IssmDouble>(size);

			ISSM_MPI_Allreduce(maxx,allmax,size,ISSM_MPI_PDOUBLE,ISSM_MPI_MAX,IssmComm::GetComm());
			ISSM_MPI_Allreduce(minx,allmin,size,ISSM_MPI_PDOUBLE,ISSM_MPI_MIN,IssmComm::GetComm());

			/*Store for later use in histogram computation:*/
			maxmeans[f]=allmax;
			minmeans[f]=allmin;

		} /*}}}*/
	}

	/*Now that we have the min and max, we can start binning. First allocate 
	 * histograms, then start filling them:*/
	IssmDouble** histogram=xNew<IssmDouble*>(nfields*nsteps);
	IssmDouble** timehistogram=xNew<IssmDouble*>(nfields);

	_printf0_("Start reading files again, this time binning values in the histogram:\n");
	/*Start opening files:*/
	for (int i=(lower_row+1);i<=upper_row;i++){
		_printf0_("reading file #: " << i << "\n");
		char file[1000];
		long int  length;
		char* buffer=NULL;

		/*string:*/
		sprintf(file,"%s/%i/%s.outbin.%i",directory,my_rank+1,model,i);

		/*open file: */
		_printf0_("    opening file:\n");
		FILE* fid=fopen(file,"rb");
		if(fid==NULL)_error_("could not open file: " << file << "\n");

		/*figure out size of file, and read the whole thing:*/
		_printf0_("    reading file:\n");
		fseek (fid, 0, SEEK_END);
		length = ftell (fid);
		fseek (fid, 0, SEEK_SET);
		buffer = xNew<char>(length);
		fread (buffer, sizeof(char), length, fid);

		/*close file:*/
		fclose (fid);

		/*create a memory stream with this buffer:*/
		_printf0_("    processing file:\n");
		fid=fmemopen(buffer, length, "rb");

		/*start reading data from the buffer directly:*/
		for (int f=0;f<nfields;f++){
			char* field=fields[f];
			fseek(fid,0,SEEK_SET);
			for (int j=0;j<nsteps;j++){
				int counter=f*nsteps+j;
				xtype[counter]=readdata(&doublemat, &doublematsize, &scalar, fid,field,steps[j]);
				if(i==(lower_row+1)){
					if(xtype[counter]==1){
						IssmDouble* localhistogram=xNewZeroInit<IssmDouble>(nbins);
						IssmDouble ma=*maxxs[counter];
						IssmDouble mi=*minxs[counter];
						int index=(scalar-mi)/(ma-mi)*nbins; if (index==nbins)index--;
						if(ma==mi)index=0;
						//_printf_( index << "|" << scalar << "|" << mi << "|" << ma << "|" << nbins << "\n");
						localhistogram[index]++;
						histogram[counter]=localhistogram;
					}
					else if (xtype[counter]==3){
						IssmDouble* localhistogram=xNewZeroInit<IssmDouble>(doublematsize*nbins);
						IssmDouble* ma=maxxs[counter];
						IssmDouble* mi=minxs[counter];
						for (int k=0;k<doublematsize;k++){
							IssmDouble scalar=doublemat[k];
							int index=(scalar-mi[k])/(ma[k]-mi[k])*nbins; if (index==nbins)index--;
							if (mi[k]==ma[k])index=0;
							_assert_(scalar<=ma[k]); _assert_(scalar>=mi[k]); _assert_(index<nbins);
							localhistogram[k*nbins+index]++;
						}
						histogram[counter]=localhistogram;
						xDelete<IssmDouble>(doublemat);
					}
					else _error_("cannot carry out statistics on type " << xtype[counter]); 
				}
				else{
					if(xtype[counter]==1){
						IssmDouble* localhistogram=histogram[counter];
						IssmDouble ma=*maxxs[counter];
						IssmDouble mi=*minxs[counter];
						int index=(scalar-mi)/(ma-mi)*nbins; if (index==nbins)index=nbins-1;
						if(ma==mi)index=0;
						localhistogram[index]++;
					}
					else if (xtype[counter]==3){
						IssmDouble* localhistogram=histogram[counter];
						IssmDouble* ma=maxxs[counter];
						IssmDouble* mi=minxs[counter];
						for (int k=0;k<doublematsize;k++){
							IssmDouble scalar=doublemat[k];
							int index=(scalar-mi[k])/(ma[k]-mi[k])*nbins; if (index==nbins)index=nbins-1;
							if (mi[k]==ma[k])index=0;

							localhistogram[k*nbins+index]++;
						}
						xDelete<IssmDouble>(doublemat);
					}
					else _error_("cannot carry out statistics on type " << xtype[counter]); 
				}
			}
		}
		_printf0_("    average in time:\n");

		/*Deal with average in time: */
		for (int f=0;f<nfields;f++){
			fseek(fid,0,SEEK_SET);
			char* field=fields[f];
			meanxtype[f]=readdata(&doublemat, &doublematsize, &scalar, fid,field,steps[0]);

			if(meanxtype[f]==1){
				IssmDouble timemean=0;
				fseek(fid,0,SEEK_SET);
				for (int j=0;j<nsteps;j++){
					readdata(&doublemat, &doublematsize, &scalar, fid,field,steps[j]);
					timemean+=scalar/nsteps;
				}

				/*Figure out max and min of time means: */
				if(i==(lower_row+1)){
					IssmDouble* localhistogram=xNewZeroInit<IssmDouble>(nbins); 
					IssmDouble ma=*maxmeans[f];
					IssmDouble mi=*minmeans[f];
					int index=(timemean-mi)/(ma-mi)*nbins; if (index==nbins)index=nbins-1;
					if(ma==mi)index=0;
					localhistogram[index]++;
					timehistogram[f]=localhistogram;
				}
				else{
					IssmDouble* localhistogram=timehistogram[f];
					IssmDouble ma=*maxmeans[f];
					IssmDouble mi=*minmeans[f];
					int index=(timemean-mi)/(ma-mi)*nbins; if (index==nbins)index=nbins-1;
					if(ma==mi)index=0;
					localhistogram[index]++;
				}
			}
			else{
				fseek(fid,0,SEEK_SET);
				IssmDouble* timemean=xNewZeroInit<IssmDouble>(doublematsize);
				for (int j=0;j<nsteps;j++){
					readdata(&doublemat, &doublematsize, &scalar, fid,field,steps[j]);
					for (int k=0;k<doublematsize;k++){
						timemean[k]+=doublemat[k]/nsteps;
					}
					xDelete<IssmDouble>(doublemat);
				}

				if(i==(lower_row+1)){
					IssmDouble* localhistogram=xNewZeroInit<IssmDouble>(doublematsize*nbins);
					IssmDouble* ma=maxmeans[f];
					IssmDouble* mi=minmeans[f];

					for (int k=0;k<doublematsize;k++){
						IssmDouble scalar=timemean[k];
						int index=(scalar-mi[k])/(ma[k]-mi[k])*nbins; if (index==nbins)index=nbins-1;
						if (mi[k]==ma[k])index=0;
						localhistogram[k*nbins+index]++;
					}
					timehistogram[f]=localhistogram;
				}
				else{

					IssmDouble* localhistogram=timehistogram[f];
					IssmDouble* ma=maxmeans[f];
					IssmDouble* mi=minmeans[f];

					for (int k=0;k<doublematsize;k++){
						IssmDouble scalar=timemean[k];
						int index=(scalar-mi[k])/(ma[k]-mi[k])*nbins; if (index==nbins)index=nbins-1;
						if (mi[k]==ma[k])index=0;

						localhistogram[k*nbins+index]++;
					}
				}
			}
		}
		fclose(fid);

		/*delete buffer:*/
		xDelete<char>(buffer);
	}
	_printf0_("Start aggregating histogram:\n");

	/*We have agregated histograms across the cluster, now gather them across  the cluster onto
	 *cpu0: */
	for (int f=0;f<nfields;f++){
		int counter0=f*nsteps+0;
		if (xtype[counter0]==1){ /*deal with scalars {{{*/
			for (int j=0;j<nsteps;j++){
				int counter=f*nsteps+j;

				/*we are broadcasting doubles:*/
				IssmDouble* histo=histogram[counter]; //size nbins
				IssmDouble* allhisto=xNewZeroInit<IssmDouble>(nbins);

				ISSM_MPI_Allreduce(histo,allhisto,nbins,ISSM_MPI_PDOUBLE,ISSM_MPI_SUM,IssmComm::GetComm());

				/*add to results:*/
				if(my_rank==0){
					char fieldname[1000];

					sprintf(fieldname,"%s%s",fields[f],"Histogram");
					results->AddResult(new GenericExternalResult<IssmPDouble*>(results->Size()+1,fieldname,allhisto,1,nbins,j+1,0));

					sprintf(fieldname,"%s%s",fields[f],"Max");
					results->AddResult(new GenericExternalResult<IssmDouble>(results->Size()+1,fieldname,*maxxs[counter],j+1,0));
					sprintf(fieldname,"%s%s",fields[f],"Min");
					results->AddResult(new GenericExternalResult<IssmDouble>(results->Size()+1,fieldname,*minxs[counter],j+1,0));
				}
			}
		} /*}}}*/
		else{ /*deal with arrays:{{{*/

			int size=xsize[counter0];
			for (int j=0;j<nsteps;j++){
				int counter=f*nsteps+j;

				/*we are broadcasting double arrays:*/
				IssmDouble* histo=histogram[counter];
				IssmDouble* allhisto=xNew<IssmDouble>(size*nbins);

				ISSM_MPI_Allreduce(histo,allhisto,size*nbins,ISSM_MPI_PDOUBLE,ISSM_MPI_SUM,IssmComm::GetComm());
				xDelete<IssmDouble>(histo);

				/*add to results:*/
				if(my_rank==0){
					char fieldname[1000];

					sprintf(fieldname,"%s%s",fields[f],"Histogram");
					results->AddResult(new GenericExternalResult<IssmPDouble*>(results->Size()+1,fieldname,allhisto,size,nbins,j+1,0));

					sprintf(fieldname,"%s%s",fields[f],"Max");
					results->AddResult(new GenericExternalResult<IssmPDouble*>(results->Size()+1,fieldname,maxxs[counter],size,1,j+1,0));
					sprintf(fieldname,"%s%s",fields[f],"Min");
					results->AddResult(new GenericExternalResult<IssmPDouble*>(results->Size()+1,fieldname,minxs[counter],size,1,j+1,0));
				}
			}
		} /*}}}*/
	}
	_printf0_("Start aggregating time mean histogram:\n");

	/*Now do the same for the time mean fields:*/
	for (int f=0;f<nfields;f++){
		if (meanxtype[f]==1){ /*deal with scalars {{{*/

			/*we are broadcasting doubles:*/
			IssmDouble* histo=timehistogram[f];
			IssmDouble* allhisto=xNewZeroInit<IssmDouble>(nbins);

			ISSM_MPI_Allreduce(histo,allhisto,nbins,ISSM_MPI_PDOUBLE,ISSM_MPI_MAX,IssmComm::GetComm());

			/*add to results at time step 1:*/
			if(my_rank==0){
				char fieldname[1000];

				sprintf(fieldname,"%s%s",fields[f],"TimeMeanHistogram");
				results->AddResult(new GenericExternalResult<IssmPDouble*>(results->Size()+1,fieldname,allhisto,1,nbins,1,0));

				sprintf(fieldname,"%s%s",fields[f],"TimeMeanMax");
				results->AddResult(new GenericExternalResult<IssmDouble>(results->Size()+1,fieldname,*maxmeans[f],1,0));
				sprintf(fieldname,"%s%s",fields[f],"TimeMeaMin");
				results->AddResult(new GenericExternalResult<IssmDouble>(results->Size()+1,fieldname,*minmeans[f],1,0));
			}
		} /*}}}*/
		else{ /*deal with arrays:{{{*/

			int size=meanxsize[f];

			/*we are broadcasting double arrays:*/
			IssmDouble* histo=timehistogram[f];
			IssmDouble* allhisto=xNewZeroInit<IssmDouble>(size*nbins);

			ISSM_MPI_Allreduce(histo,allhisto,size*nbins,ISSM_MPI_PDOUBLE,ISSM_MPI_SUM,IssmComm::GetComm());
			xDelete<IssmDouble>(histo);
			/*add to results at step 1:*/
			if(my_rank==0){
				char fieldname[1000];

				sprintf(fieldname,"%s%s",fields[f],"TimeMeanHistogram");
				results->AddResult(new GenericExternalResult<IssmPDouble*>(results->Size()+1,fieldname,allhisto,size,nbins,1,0));
				sprintf(fieldname,"%s%s",fields[f],"TimeMeanMax");
				results->AddResult(new GenericExternalResult<IssmPDouble*>(results->Size()+1,fieldname,maxmeans[f],size,1,1,0));
				sprintf(fieldname,"%s%s",fields[f],"TimeMeanMin");
				results->AddResult(new GenericExternalResult<IssmPDouble*>(results->Size()+1,fieldname,minmeans[f],size,1,1,0));
			}
		} /*}}}*/
	}
	_printf0_("Done aggregating time mean histogram:\n");
	IssmComm::SetComm(ISSM_MPI_COMM_WORLD);

	return 1;
}
/*}}}*/
int ComputeMeanVariance(Parameters* parameters,Results* results,int color, ISSM_MPI_Comm statcomm){  /*{{{*/

	int nsamples; 
	char* directory=NULL;
	char* model=NULL;
	char** fields=NULL;
	int* steps=NULL;
	int nsteps;
	int nfields;
	int range,lower_row,upper_row;
	int nfilesperdirectory;

	/*intermediary:*/
	IssmDouble* doublemat=NULL;
	int         doublematsize;
	IssmDouble scalar;

	/*computation of average and variance itself:*/
	IssmDouble*  x = NULL;
	IssmDouble*  x2 = NULL;
	IssmDouble** xs = NULL;
	IssmDouble** xs2 = NULL;
	int*         xtype=NULL;
	int*         xsize=NULL;

	IssmDouble** meanx=NULL;
	IssmDouble** meanx2=NULL;
	int*         meantype=NULL;
	int*         meansize=NULL;

	/*only work on the statistical communicator: */
	if (color==MPI_UNDEFINED)return 0;

	/*Retrieve parameters:*/
	parameters->FindParam(&nfilesperdirectory,QmuNfilesPerDirectoryEnum);
	parameters->FindParam(&nsamples,QmuNsampleEnum);
	parameters->FindParam(&directory,DirectoryNameEnum);
	parameters->FindParam(&model,InputFileNameEnum);
	parameters->FindParam(&fields,&nfields,FieldsEnum);
	parameters->FindParam(&steps,&nsteps,StepsEnum);

	/*Get rank from the stat comm communicator:*/
	IssmComm::SetComm(statcomm);
	int my_rank=IssmComm::GetRank();

	/*Open files and read them complelety, in a distributed way:*/
	range=DetermineLocalSize(nsamples,IssmComm::GetComm());
	GetOwnershipBoundariesFromRange(&lower_row,&upper_row,range,IssmComm::GetComm());

	/*Initialize arrays:*/
	xs=xNew<IssmDouble*>(nfields*nsteps);
	xs2=xNew<IssmDouble*>(nfields*nsteps);
	xtype=xNew<int>(nfields*nsteps);
	xsize=xNew<int>(nfields*nsteps);

	meantype=xNew<int>(nfields);
	meansize=xNew<int>(nfields);
	meanx=xNew<IssmDouble*>(nfields);
	meanx2=xNew<IssmDouble*>(nfields);

	/*Start opening files:*/
	for (int i=(lower_row+1);i<=upper_row;i++){
		_printf0_("reading file #: " << i << "\n");
		char file[1000];
		long int  length;
		char* buffer=NULL;

		/*string:*/
		sprintf(file,"%s/%i/%s.outbin.%i",directory,my_rank+1,model,i);

		/*open file: */
		_printf0_("    opening file: " << file << "\n");
		FILE* fid=fopen(file,"rb");
		if(fid==NULL) _error_("    could not open file: " << file << "\n");

		/*figure out size of file, and read the whole thing:*/
		_printf0_("    reading file:\n");
		fseek (fid, 0, SEEK_END);
		length = ftell (fid);
		fseek (fid, 0, SEEK_SET);
		buffer = xNew<char>(length);
		fread (buffer, sizeof(char), length, fid);

		/*close file:*/
		fclose (fid);

		/*create a memory stream with this buffer:*/
		_printf0_("    processing file:\n");
		fid=fmemopen(buffer, length, "rb");

		/*start reading data from the buffer directly:*/
		for (int f=0;f<nfields;f++){
			char* field=fields[f];
			fseek(fid,0,SEEK_SET);
			for (int j=0;j<nsteps;j++){
				int counter=f*nsteps+j;
				xtype[counter]=readdata(&doublemat, &doublematsize, &scalar, fid,field,steps[j]);
				if(i==(lower_row+1)){
					if(xtype[counter]==1){
						xs[counter]=xNew<IssmDouble>(1); 
						xs2[counter]=xNew<IssmDouble>(1); 
						*xs[counter]=scalar;
						*xs2[counter]=pow(scalar,2.0);
						xsize[counter]=1;
					}
					else if (xtype[counter]==3){
						IssmDouble* doublemat2=xNew<IssmDouble>(doublematsize);
						for(int k=0;k<doublematsize;k++)doublemat2[k]=pow(doublemat[k],2.0);
						xs[counter]=doublemat;
						xs2[counter]=doublemat2;
						xsize[counter]=doublematsize;
					}
					else _error_("cannot carry out statistics on type " << xtype[counter]); 
				}
				else{
					if(xtype[counter]==1){
						*xs[counter]+=scalar;
						*xs2[counter]+=pow(scalar,2.0);
					}
					else if (xtype[counter]==3){
						IssmDouble* newdoublemat=xs[counter];
						IssmDouble* newdoublemat2=xs2[counter];
						for(int k=0;k<doublematsize;k++){
							newdoublemat[k]+=doublemat[k];
							newdoublemat2[k]+=pow(doublemat[k],2.0);
						}
						xs[counter]=newdoublemat;
						xs2[counter]=newdoublemat2;
					}
					else _error_("cannot carry out statistics on type " << xtype[counter]); 
				}
			}
		}

		/*Deal with time mean: */
		for (int f=0;f<nfields;f++){
			char* field=fields[f];
			fseek(fid,0,SEEK_SET);
			meantype[f]=readdata(&doublemat, &doublematsize, &scalar, fid,field,steps[0]);
			if(i==(lower_row+1)){
				if(meantype[f]==1){
					meanx[f]=xNewZeroInit<IssmDouble>(1);
					meanx2[f]=xNewZeroInit<IssmDouble>(1);
					meansize[f]=1;
				}
				else{
					meanx[f]=xNewZeroInit<IssmDouble>(doublematsize);
					meanx2[f]=xNewZeroInit<IssmDouble>(doublematsize);
					meansize[f]=doublematsize;
				}
			}
			fseek(fid,0,SEEK_SET);
			if(meantype[f]==1){
				IssmDouble sc=0;
				IssmDouble sc2=0;
				for(int j=0;j<nsteps;j++){
					readdata(&doublemat, &doublematsize, &scalar, fid,field,steps[j]);
					sc+=scalar/nsteps;
				}
				sc2+=pow(sc,2.0);
				*meanx[f]+=sc;
				*meanx2[f]+=sc2;
			}
			else{
				IssmDouble* sc=meanx[f];
				IssmDouble* sc2=meanx2[f];
				IssmDouble* timemean=xNewZeroInit<IssmDouble>(doublematsize);
				IssmDouble* timemean2=xNewZeroInit<IssmDouble>(doublematsize);

				for(int j=0;j<nsteps;j++){
					readdata(&doublemat, &doublematsize, &scalar, fid,field,steps[j]);
					for (int k=0;k<doublematsize;k++){
						timemean[k]+=doublemat[k]/nsteps;
					}
				}
				for (int k=0;k<doublematsize;k++){
					timemean2[k]=pow(timemean[k],2.0);
				}
				for (int k=0;k<doublematsize;k++){
					sc[k]+=timemean[k];
					sc2[k]+=timemean2[k];
				}

			}

		}
		fclose(fid);

		/*delete buffer:*/
		xDelete<char>(buffer);
	}
	ISSM_MPI_Barrier(IssmComm::GetComm());
	_printf0_("Done reading files, now computing mean and variance.\n"); 

	/*We have agregated x and x^2 across the cluster, now gather across the cluster onto
	 *cpu0 and then compute statistics:*/
	for (int f=0;f<nfields;f++){
		int counter0=f*nsteps+0;
		if (xtype[counter0]==1){ /*deal with scalars {{{*/
			IssmDouble mean,stddev;
			for (int j=0;j<nsteps;j++){
				int counter=f*nsteps+j;

				/*we are broadcasting doubles:*/
				IssmDouble scalar=*xs[counter];
				IssmDouble scalar2=*xs2[counter];
				IssmDouble sumscalar;
				IssmDouble sumscalar2;

				ISSM_MPI_Reduce(&scalar,&sumscalar,1,ISSM_MPI_PDOUBLE,ISSM_MPI_SUM,0,IssmComm::GetComm());
				ISSM_MPI_Reduce(&scalar2,&sumscalar2,1,ISSM_MPI_PDOUBLE,ISSM_MPI_SUM,0,IssmComm::GetComm());
				/*Build average and standard deviation. For standard deviation, use the 
				 *following formula: sigma^2=E(x^2)-mu^2:*/
				mean=sumscalar/nsamples;
				stddev=sqrt(sumscalar2/nsamples-pow(mean,2.0));

				/*add to results:*/
				if(my_rank==0){
					char fieldname[1000];

					sprintf(fieldname,"%s%s",fields[f],"Mean");
					results->AddResult(new GenericExternalResult<IssmDouble>(results->Size()+1,fieldname,mean,j+1,0));
					sprintf(fieldname,"%s%s",fields[f],"Stddev");
					results->AddResult(new GenericExternalResult<IssmDouble>(results->Size()+1,fieldname,stddev,j+1,0));
				}

			}
		} /*}}}*/
		else{ /*deal with arrays:{{{*/

			int size=xsize[counter0];

			IssmDouble*  mean=xNew<IssmDouble>(size);
			IssmDouble*  stddev=xNew<IssmDouble>(size);

			for (int j=0;j<nsteps;j++){
				int counter=f*nsteps+j;

				/*we are broadcasting double arrays:*/
				x=xs[counter];
				x2=xs2[counter];

				IssmDouble*  sumx=xNew<IssmDouble>(size);
				IssmDouble*  sumx2=xNew<IssmDouble>(size);

				ISSM_MPI_Reduce(x,sumx,size,ISSM_MPI_PDOUBLE,ISSM_MPI_SUM,0,IssmComm::GetComm());
				ISSM_MPI_Reduce(x2,sumx2,size,ISSM_MPI_PDOUBLE,ISSM_MPI_SUM,0,IssmComm::GetComm());

				/*Build average and standard deviation. For standard deviation, use the 
				 *following formula: sigma^2=E(x^2)-mu^2:*/
				for (int k=0;k<size;k++){
					mean[k]=sumx[k]/nsamples;
					stddev[k]=sqrt(sumx2[k]/nsamples-pow(mean[k],2.0));
				}

				/*add to results:*/
				if(my_rank==0){
					char fieldname[1000];

					sprintf(fieldname,"%s%s",fields[f],"Mean");
					results->AddResult(new GenericExternalResult<IssmPDouble*>(results->Size()+1,fieldname,mean,size,1,j+1,0));
					sprintf(fieldname,"%s%s",fields[f],"Stddev");
					results->AddResult(new GenericExternalResult<IssmPDouble*>(results->Size()+1,fieldname,stddev,size,1,j+1,0));
				}
			}
		} /*}}}*/
	}
	/*Do the same but for the time mean:*/
	for (int f=0;f<nfields;f++){
		if (meantype[f]==1){ /*deal with scalars {{{*/
			IssmDouble mean,stddev;

			/*we are broadcasting doubles:*/
			IssmDouble scalar=*meanx[f];
			IssmDouble scalar2=*meanx2[f];
			IssmDouble sumscalar;
			IssmDouble sumscalar2;

			ISSM_MPI_Reduce(&scalar,&sumscalar,1,ISSM_MPI_PDOUBLE,ISSM_MPI_SUM,0,IssmComm::GetComm());
			ISSM_MPI_Reduce(&scalar2,&sumscalar2,1,ISSM_MPI_PDOUBLE,ISSM_MPI_SUM,0,IssmComm::GetComm());
			/*Build average and standard deviation. For standard deviation, use the 
			 *following formula: sigma^2=E(x^2)-mu^2:*/
			mean=sumscalar/nsamples;
			stddev=sqrt(sumscalar2/nsamples-pow(mean,2.0));

			/*add to results:*/
			if(my_rank==0){
				char fieldname[1000];

				sprintf(fieldname,"%s%s",fields[f],"TimeMean");
				results->AddResult(new GenericExternalResult<IssmDouble>(results->Size()+1,fieldname,mean,1,0));
				sprintf(fieldname,"%s%s",fields[f],"TimeStddev");
				results->AddResult(new GenericExternalResult<IssmDouble>(results->Size()+1,fieldname,stddev,1,0));
			}
		} /*}}}*/
		else{ /*deal with arrays:{{{*/

			int size=meansize[f];
			IssmDouble*  mean=xNew<IssmDouble>(size);
			IssmDouble*  stddev=xNew<IssmDouble>(size);

			/*we are broadcasting double arrays:*/
			x=meanx[f];
			x2=meanx2[f];

			IssmDouble*  sumx=xNew<IssmDouble>(size);
			IssmDouble*  sumx2=xNew<IssmDouble>(size);

			ISSM_MPI_Reduce(x,sumx,size,ISSM_MPI_PDOUBLE,ISSM_MPI_SUM,0,IssmComm::GetComm());
			ISSM_MPI_Reduce(x2,sumx2,size,ISSM_MPI_PDOUBLE,ISSM_MPI_SUM,0,IssmComm::GetComm());

			/*Build average and standard deviation. For standard deviation, use the 
			 *following formula: sigma^2=E(x^2)-mu^2:*/
			for (int k=0;k<size;k++){
				mean[k]=sumx[k]/nsamples;
				stddev[k]=sqrt(sumx2[k]/nsamples-pow(mean[k],2.0));
			}

			/*add to results:*/
			if(my_rank==0){
				char fieldname[1000];

				sprintf(fieldname,"%s%s",fields[f],"TimeMean");
				results->AddResult(new GenericExternalResult<IssmPDouble*>(results->Size()+1,fieldname,mean,size,1,1,0));
				sprintf(fieldname,"%s%s",fields[f],"TimeStddev");
				results->AddResult(new GenericExternalResult<IssmPDouble*>(results->Size()+1,fieldname,stddev,size,1,1,0));
			}
		} /*}}}*/
	}

	_printf0_("Done with MeanVariance:\n");
	IssmComm::SetComm(ISSM_MPI_COMM_WORLD);

	return 1;
} /*}}}*/
int ComputeSampleSeries(Parameters* parameters,Results* results,int color, ISSM_MPI_Comm statcomm){ /*{{{*/

	int nsamples; 
	char* directory=NULL;
	char* model=NULL;
	char** fields=NULL;
	int* steps=NULL;
	int nsteps;
	int nfields;
	int range,lower_row,upper_row;
	int nfilesperdirectory;
	int* indices=NULL;
	int  nindices;

	/*intermediary:*/
	IssmDouble* doublemat=NULL;
	int         doublematsize;
	IssmDouble scalar;

	/*computation of average and variance itself:*/
	IssmDouble*  x = NULL;
	IssmDouble*  allx=NULL;
	IssmDouble** xs = NULL;
	int*         xtype=NULL;
	int*         xsize=NULL;

	/*only work on the statistical communicator: */
	if (color==MPI_UNDEFINED)return 0;

	/*Retrieve parameters:*/
	parameters->FindParam(&nsamples,QmuNsampleEnum);
	parameters->FindParam(&directory,DirectoryNameEnum);
	parameters->FindParam(&model,InputFileNameEnum);
	parameters->FindParam(&fields,&nfields,FieldsEnum);
	parameters->FindParam(&steps,&nsteps,StepsEnum);
	parameters->FindParam(&indices,&nindices,IndicesEnum);

	/*Get rank from the stat comm communicator:*/
	IssmComm::SetComm(statcomm);
	int my_rank=IssmComm::GetRank();

	/*Open files and read them complelety, in a distributed way:*/
	range=DetermineLocalSize(nsamples,IssmComm::GetComm());
	GetOwnershipBoundariesFromRange(&lower_row,&upper_row,range,IssmComm::GetComm());

	/*Initialize arrays:*/
	xs=xNew<IssmDouble*>(nfields*nsteps);
	xtype=xNew<int>(nfields*nsteps);
	xsize=xNew<int>(nfields*nsteps);

	/*Start opening files:*/
	for (int i=(lower_row+1);i<=upper_row;i++){
		_printf0_("reading file #: " << i << "\n");
		char file[1000];
		long int  length;
		char* buffer=NULL;

		/*string:*/
		sprintf(file,"%s/%i/%s.outbin.%i",directory,my_rank+1,model,i);

		/*open file: */
		_printf0_("    opening file:\n");
		FILE* fid=fopen(file,"rb");

		/*figure out size of file, and read the whole thing:*/
		_printf0_("    reading file:\n");
		fseek (fid, 0, SEEK_END);
		length = ftell (fid);
		fseek (fid, 0, SEEK_SET);
		buffer = xNew<char>(length);
		fread (buffer, sizeof(char), length, fid);

		/*close file:*/
		fclose (fid);

		/*create a memory stream with this buffer:*/
		_printf0_("    processing file:\n");
		fid=fmemopen(buffer, length, "rb");

		/*start reading data from the buffer directly:*/
		for (int f=0;f<nfields;f++){
			fseek(fid,0,SEEK_SET);
			char* field=fields[f];
			for (int j=0;j<nsteps;j++){
				int counter=f*nsteps+j;
				xtype[counter]=readdata(&doublemat, &doublematsize, &scalar, fid,field,steps[j]);
				if(i==(lower_row+1)){
					if(xtype[counter]==1){
						x=xNew<IssmDouble>(range);
						x[0]=scalar;
						xs[counter]=x;
						xsize[counter]=range;
					}
					else if (xtype[counter]==3){
						x=xNew<IssmDouble>(nindices*range);
						for(int k=0;k<nindices;k++)x[(i-(lower_row+1))*nindices+k]=doublemat[indices[k]-1];
						xs[counter]=x;
						xsize[counter]=range*nindices;
					}
					else _error_("cannot carry out statistics on type " << xtype[counter]); 
				}
				else{
					if(xtype[counter]==1){
						x=xs[counter]; 
						x[i-(lower_row+1)]=scalar;
						xs[counter]=x;
					}
					else if (xtype[counter]==3){
						x=xs[counter];
						for(int k=0;k<nindices;k++)x[(i-(lower_row+1))*nindices+k]=doublemat[indices[k]-1];
						xs[counter]=x;
					}
					else _error_("cannot carry out statistics on type " << xtype[counter]); 
				}
			}
		}
		fclose(fid);

		/*delete buffer:*/
		xDelete<char>(buffer);
	}
	ISSM_MPI_Barrier(IssmComm::GetComm());
	_printf0_("Done reading files, now assembling time series.\n");

	for (int f=0;f<nfields;f++){
		for (int j=0;j<nsteps;j++){
			int counter=f*nsteps+j;
			if (xtype[counter]==1){
				/*we are broadcasting range times doubles:*/
				x=xs[counter]; 
				allx=xNew<IssmDouble>(nsamples);
				MPI_Gather(x, range, ISSM_MPI_PDOUBLE,allx, range, ISSM_MPI_PDOUBLE, 0, IssmComm::GetComm());
				/*add to results:*/
				if(my_rank==0){
					char fieldname[1000];

					sprintf(fieldname,"%s%s",fields[f],"Samples");
					results->AddResult(new GenericExternalResult<IssmPDouble*>(results->Size()+1,fieldname,allx,nsamples,1,j+1,0));
				}
			}
			else{
				/*we are broadcasting double arrays:*/
				x=xs[counter];
				allx=xNew<IssmDouble>(nsamples*nindices);

				MPI_Gather(x, range*nindices, ISSM_MPI_PDOUBLE,allx, range*nindices, ISSM_MPI_PDOUBLE, 0, IssmComm::GetComm());

				/*add to results:*/
				if(my_rank==0){
					char fieldname[1000];
					sprintf(fieldname,"%s%s",fields[f],"Samples");
					results->AddResult(new GenericExternalResult<IssmPDouble*>(results->Size()+1,fieldname,allx,nsamples,nindices,j+1,0));
				}
			}
		}
	}
	_printf0_("Done with SampleSeries:\n");
	IssmComm::SetComm(ISSM_MPI_COMM_WORLD);

	return 1;
} /*}}}*/
int OutputStatistics(Parameters* parameters,Results* results,int color,ISSM_MPI_Comm statcomm){ /*{{{*/

	char   outputfilename[1000];
	char* directory=NULL;
	char* model=NULL;
	char* method=NULL;
	int   nsamples;
	int* steps=NULL;
	int nsteps;

	/*only work on the statistical communicator: */
	if (color==MPI_UNDEFINED)return 0;

	FemModel* femmodel=new FemModel();

	/*Some parameters that will allow us to use the OutputResultsx module:*/
	parameters->AddObject(new BoolParam(QmuIsdakotaEnum,false));
	parameters->AddObject(new BoolParam(SettingsIoGatherEnum,true));

	parameters->FindParam(&directory,DirectoryNameEnum);
	parameters->FindParam(&model,InputFileNameEnum);
	parameters->FindParam(&nsamples,QmuNsampleEnum);
	parameters->FindParam(&steps,&nsteps,StepsEnum);

	sprintf(outputfilename,"%s/%s.stats",directory,model);
	parameters->AddObject(new StringParam(OutputFileNameEnum,outputfilename));

	/*Call OutputResults module:*/
	femmodel->parameters=parameters;
	femmodel->results=results;

	OutputResultsx(femmodel);

	return 1;
} /*}}}*/
bool DakotaDirStructure(int argc,char** argv){ /*{{{*/

	char* input_file; 
	FILE* fid;
	IoModel* iomodel=NULL;
	int check;

	//qmu statistics
	bool statistics    = false;
	int  numdirectories = 0;

	/*First things first, set the communicator as a global variable: */
	IssmComm::SetComm(MPI_COMM_WORLD);

	/*Barrier:*/
	ISSM_MPI_Barrier(IssmComm::GetComm());
	_printf0_("Preparing directory structure for model outputs:" << "\n");

	//open model input file for reading
	input_file=xNew<char>((strlen(argv[2])+strlen(argv[3])+strlen(".bin")+2));
	sprintf(input_file,"%s/%s%s",argv[2],argv[3],".bin");
	fid=fopen(input_file,"rb");
	if (fid==NULL) Cerr << "dirstructure error message: could not open model " << input_file << " to retrieve qmu statistics parameters" << std::endl;

	//initialize IoModel, but light version, we just need it to fetch one constant: 
	iomodel=new IoModel();
	iomodel->fid=fid;
	iomodel->FetchConstants();

	//early return if statistics not requested: 
	iomodel->FindConstant(&statistics,"md.qmu.statistics");
	if(!statistics){
		delete iomodel;
		xDelete<char>(input_file);
		fclose(fid);
		return false; //important return value!
	}

	iomodel->FindConstant(&numdirectories,"md.qmu.statistics.ndirectories");

	/*Ok, we have everything we need to create the directory structure:*/
	if(IssmComm::GetRank()==0){
		for (int i=0;i<numdirectories;i++){
			char directory[1000];
			sprintf(directory,"./%i",i+1);

			check = mkdir(directory,ACCESSPERMS);
			if (check) _error_("dirstructure error message: could not create directory " << directory << "\n");
		}
	}

	/*Delete resources:*/
	delete iomodel;
	xDelete<char>(input_file);

	//close model file: 
	fclose(fid);

	//return value: 
	return true; //statistics computation on!
} /*}}}*/
int DakotaStatistics(int argc,char** argv){ /*{{{*/

	char* input_file; 
	FILE* fid;
	IoModel* iomodel=NULL;
	ISSM_MPI_Comm statcomm;
	int my_rank;

	//qmu statistics
	bool statistics    = false;
	int  numstatistics = 0;
	int  numdirectories = 0;
	int  nfilesperdirectory = 0;
	char string[1000];
	char* name = NULL;
	char** fields = NULL;
	int    nfields; 
	int*   steps=NULL;
	int    nsteps;
	int    nbins;
	int*   indices=NULL;
	int    nindices;
	int    nsamples;
	int    dummy;
	char*  directory=NULL;
	char*  model=NULL;
	Results* results=NULL;
	Parameters* parameters=NULL;
	int color;

	/*First things first, set the communicator as a global variable: */
	IssmComm::SetComm(MPI_COMM_WORLD);
	my_rank=IssmComm::GetRank();

	/*Barrier:*/
	ISSM_MPI_Barrier(IssmComm::GetComm());
	_printf0_("Dakota Statistic Computation" << "\n");

	//open model input file for reading
	input_file=xNew<char>((strlen(argv[2])+strlen(argv[3])+strlen(".bin")+2));
	sprintf(input_file,"%s/%s%s",argv[2],argv[3],".bin");
	fid=fopen(input_file,"rb");
	if (fid==NULL) Cerr << "issm_dakota_statistics error message: could not open model " << input_file << " to retrieve qmu statistics parameters" << std::endl;

	//initialize IoModel, but light version, we'll need it to fetch constants:
	iomodel=new IoModel();
	iomodel->fid=fid;
	iomodel->FetchConstants();

	//early return if statistics not requested: 
	iomodel->FindConstant(&statistics,"md.qmu.statistics");
	if(!statistics){
		delete iomodel;
		xDelete<char>(input_file);
		fclose(fid); 
		return 0;
	}else{
		//create parameters datasets with al the qmu statistics settings we need: 

		/*Initialize parameters and results:*/
		results   = new Results();
		parameters=new Parameters();

		//solution type: 
		parameters->AddObject(new IntParam(SolutionTypeEnum,StatisticsSolutionEnum));

		//root  directory
		directory=xNew<char>(strlen(argv[2])+1);
		xMemCpy<char>(directory,argv[2],strlen(argv[2])+1);
		parameters->AddObject(new StringParam(DirectoryNameEnum,directory));

		//model  name
		model=xNew<char>(strlen(argv[3])+1);
		xMemCpy<char>(model,argv[3],strlen(argv[3])+1);
		parameters->AddObject(new StringParam(InputFileNameEnum,model));

		//nsamples
		iomodel->FindConstant(&nsamples,"md.qmu.method.params.samples");
		parameters->AddObject(new IntParam(QmuNsampleEnum,nsamples));

		//ndirectories
		iomodel->FindConstant(&numdirectories,"md.qmu.statistics.ndirectories");
		parameters->AddObject(new IntParam(QmuNdirectoriesEnum,numdirectories));

		//nfiles per directory
		iomodel->FindConstant(&nfilesperdirectory,"md.qmu.statistics.nfiles_per_directory");
		parameters->AddObject(new IntParam(QmuNfilesPerDirectoryEnum,nfilesperdirectory));

		//At this point, we don't want to go forward any longer, we want to create an MPI 
		//communicator on which to carry out the computations:
		if ((my_rank+1)*nfilesperdirectory>nsamples)color=MPI_UNDEFINED;
		else color=0;
		ISSM_MPI_Comm_split(ISSM_MPI_COMM_WORLD,color, my_rank, &statcomm);

		iomodel->FindConstant(&numstatistics,"md.qmu.statistics.numstatistics");
		for (int i=1;i<=numstatistics;i++){

			char* directory=NULL;
			char* model=NULL;
			int   nsamples;
			_printf0_("Dealing with qmu statistical computation #" << i << "\n");

			sprintf(string,"md.qmu.statistics.method(%i).name",i);
			iomodel->FindConstant(&name,string);

			sprintf(string,"md.qmu.statistics.method(%i).fields",i);
			iomodel->FindConstant(&fields,&nfields,string);
			parameters->AddObject(new StringArrayParam(FieldsEnum,fields,nfields));

			sprintf(string,"md.qmu.statistics.method(%i).steps",i);
			iomodel->FetchData(&steps,&dummy,&nsteps,string);
			parameters->AddObject(new IntVecParam(StepsEnum,steps,nsteps));

			if (strcmp(name,"Histogram")==0){
				/*fetch nbins: */
				sprintf(string,"md.qmu.statistics.method(%i).nbins",i);
				iomodel->FindConstant(&nbins,string);
				parameters->AddObject(new IntParam(NbinsEnum,nbins));
				ComputeHistogram(parameters,results,color,statcomm);
			}
			else if (strcmp(name,"SampleSeries")==0){
				/*fetch indices: */
				sprintf(string,"md.qmu.statistics.method(%i).indices",i);
				iomodel->FetchData(&indices,&dummy,&nindices,string);
				parameters->AddObject(new IntVecParam(IndicesEnum,indices,nindices));

				ComputeSampleSeries(parameters,results,color,statcomm);
			}
			else if (strcmp(name,"MeanVariance")==0){
				ComputeMeanVariance(parameters,results,color,statcomm);
			}
			else _error_(" error creating qmu statistics methods parameters: unsupported method " << name);
		}

		/*Delete resources:*/
		xDelete<char>(input_file);
		delete iomodel;
	}

	//close model file: 
	fclose(fid);

	/*output results:*/
	OutputStatistics(parameters,results,color,statcomm);

	/*all meet here: */
	ISSM_MPI_Barrier(ISSM_MPI_COMM_WORLD); _printf0_("Output file.\n");

	/*Delete resources:*/
	delete parameters; 
	delete results;

	return 1;
} /*}}}*/
