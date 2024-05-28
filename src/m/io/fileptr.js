//FILEPTR class definition
////
//// Usage: 
//// var fid = new fileptr(); 
//
//

function fileptr() {
	//properties
	this.increment=NaN;
	this.buffer =NaN;
	this.view =NaN;
	this.ptr =NaN;
	this.buffersize =NaN;
	this.mode ='';
	this.options = new pairoptions(Array.prototype.slice.call(arguments));
	
	//methods
		this.disp = function () { //{{{
			console.log(sprintf("   fileptr:")); 

			console.log(sprintf("       buffer: ArrayBuffer{ byteLength: %i }\n",this.buffer.byteLength));
			console.log(sprintf("       ptr: %i\n",this.ptr));
			console.log(sprintf("       increment: %i\n",this.increment));
			console.log(sprintf("       mode: %s\n",this.mode));

		} //}}}
		this.setdefaultparameters = function (options) { //{{{
	
			this.mode=options.getfieldvalue('mode');
			this.ptr=0;
			this.increment=0;
			this.buffersize=0;
			if (this.mode=='w'){
				this.increment=options.getfieldvalue('increment',8000000); //80000 bytes,  10000 doubles.
				this.buffer=new ArrayBuffer(this.increment);
				this.view=new DataView(this.buffer);
			}
			else if(this.mode == 'r'){
				
				/*recover buffer and its size: */
				var bufferin= options.getfieldvalue('buffer');
				this.buffersize= options.getfieldvalue('buffersize');
				
				/*crete a typed array buffer: */
				this.buffer=new ArrayBuffer(this.buffersize);
				this.view=new DataView(this.buffer); 
				for(var i=0;i<this.buffersize;i++) this.view.setUint8(i,bufferin[i]);
			}

		} //}}}
		this.fwrite = function (value,format) { //{{{

			
			if(format == 'int'){
				if(this.ptr+4>=this.buffer.byteLength)this.resize();
				this.view.setUint32(this.ptr,value,true); this.ptr+=4;
			}
			else if(format == 'char'){
				if(this.ptr+value.length>=this.buffer.byteLength)this.resize();
				for(var i=0;i<value.length;i++){
					this.view.setUint8(this.ptr,value.charCodeAt(i),true); 
					this.ptr+=1;
				}
			}
			else if(format == 'double'){
				if(this.ptr+8>=this.buffer.byteLength)this.resize();
				if (!IsArray(value)){
					this.view.setFloat64(this.ptr,value,true);
					this.ptr+=8;
				}
				else{
					if (!IsArray(value[0])){
						if(this.ptr+value.length*8>=this.buffer.byteLength){
							this.resize();
							if(this.ptr+value.length*8>=this.buffer.byteLength)throw Error('fileptr.fwrite error: need to increase increment size!');
						}
						for(var i=0;i<value.length;i++){
							this.view.setFloat64(this.ptr,value[i],true);
							this.ptr+=8;
						}
					}
					else{
						if(this.ptr+value.length*value[0].length*8>=this.buffer.byteLength)this.resize();
						for(var i=0;i<value.length;i++){
							for(var j=0;j<value[0].length;j++){
								this.view.setFloat64(this.ptr,value[i][j],true);
								this.ptr+=8;
							}
						}
					}
				}
			}
			else throw Error('fileptr.fwrite error message: wrong type of format');
		} //}}}
		this.fread = function (size,format) { //{{{
			
			var value;

			if(this.ptr==this.buffersize)return -1;
			if(format == 'int'){
				if(size==1){
					value=this.view.getInt32(this.ptr,true); 
					this.ptr+=4;
				}
				else{
					value = new Int32Array(size);
					for(var i=0;i<size;i++){
						value[i]=this.view.getInt32(this.ptr,true); 
						this.ptr+=4;
					}
				}
			}
			else if(format == 'char'){
				value = ''; 
				for(var i=0;i<(size-1);i++){
					value+= String.fromCharCode(this.view.getUint8(this.ptr,true));
					this.ptr+=1;
				}
				this.ptr+=1; //pass over the '\0';

			}
			else if(format == 'double'){
				if(size==1){
					value=this.view.getFloat64(this.ptr,true);
					this.ptr+=8;
				}
				else{ 
					value = new Float64Array(size);
					for(var i=0;i<size;i++){
						value[i]=this.view.getFloat64(this.ptr,true);
						this.ptr+=8;
					}
				}
			}
			else throw Error('fileptr.fwrite error message: wrong type of format');
			
			return value;
		} //}}}
		this.rawbuffer = function () { //{{{
			return this.buffer.slice(0,this.ptr);
		} //}}}
		this.resize = function () { //{{{
			var  newbuffer = new ArrayBuffer(this.buffer.byteLength+this.increment);
			new Uint8Array(newbuffer).set(new Uint8Array(this.buffer));
			this.buffer=newbuffer;
			this.view=new DataView(this.buffer);
		} //}}}
	//set defaults
	this.setdefaultparameters(this.options);
}
