function sliderInit(){ //{{{
	/*
	Name:
		sliderInit
	
	Description:
		Initialize slider corresponding to passed selector with passed value,
		minimum and maximum values, step (which is the increment for the 
		slider), and callback function (which is called when slider is changed). 
	
		Options are passed as pairs to the function according to the usage 
		example below.
	
	Usage:
		sliderInit(
			'selector', #some-element, 
			'value', 0,
			'min', -1,
			'max', 1,
			'step', 0.1,
			'callback', function(value){someEngineFunction(value)}
		);
	
	NOTE:	jQuery Mobile will reflect changes to slider in its input "label" 
			based on the width of the slider. For example, step might be set to
			0.1, but if the width of the slider is too narrow, changes to the 
			slider will only be reflected in whole numeral increments in the
			slider' label. That said, *all* sliders can be adjusted via the
			"up" and "down" keys by whatever increment step is set to.
	*/

	// Convert arguments to options
	var args 		= Array.prototype.slice.call(arguments);
	var options 	= new pairoptions(args.slice());

	// Recover option values
	var selector	= options.getfieldvalue('selector', '');
	var value 		= options.getfieldvalue('value', 0);
	var callback 	= options.getfieldvalue('callback', function(event, ui){});
	var min 		= options.getfieldvalue('min', 0.6 * value);
	var max 		= options.getfieldvalue('max', 1.4 * value);
	var step 		= options.getfieldvalue('step', 1);

	/*
		Update slider attributes.
		
		NOTE:	Although slider has already been created, need to call slider() 
				in order to avoid:
				
					Error: cannot call methods on slider prior to 
					initialization; attempted to call method 'refresh'
					
				Attempted all other methods for intialization of slider widget
				from jQuery Mobile, and this is the only one that seemed to work
				(see index.php for related markup).
	*/
	$(selector).slider();
	$(selector).val(value);
	$(selector).attr('min', min);
	$(selector).attr('max', max);
	$(selector).attr('step', step);
	$(selector).on('slidestop', function(event, ui){
		callback(parseFloat($(selector).val()));
	});
	$(selector).slider('refresh'); //Slider must be "refreshed" after any JavaScript change to it, as it is an AJAX object.
} //}}}

function sliderMoveInput(selector){ //{{{
	/*
	Name:
		sliderMoveInput
	
	Description:
		Appends a jQuery Mobile slider input to an element whose selector 
		adheres to the following protocol,
		
			destination = sliderSelector + '-value'
			
	Usage:
		sliderMoveInput('#someSliderSelector');
		
	NOTE:	Destination element must, obviously, be hardcoded into markup for a
			call to this function to work as expected.
	*/
	$(selector).appendTo(selector + '-value');
	$(selector).slider('refresh'); //Slider must be "refreshed" after any JavaScript change to it, as it is an AJAX object.
} //}}}

function progressInit(){ //{{{
	
	// Convert arguments to options.
	var args 			= Array.prototype.slice.call(arguments);
	var options 		= new pairoptions(args.slice());
	
	// Recover option values
	var sim 			= options.getfieldvalue('sim', '');
	
	var canvas 			= $(sim + '-canvas')[0];
	var progressBar 	= $(sim + '-controls-slider-progress');
	var playButton 		= $(sim + '-controls-button-play');
	var reverseButton 	= $(sim + '-controls-button-reverse');
	var timeText 		= $(sim + '-controls-progress-time');
	
	/*
		Update slider attributes.
		
		NOTE:	Although slider has already been created, need to call slider() 
				in order to avoid:
				
					Error: cannot call methods on slider prior to 
					initialization; attempted to call method 'refresh'
					
				Attempted all other methods for intialization of slider widget
				from jQuery Mobile, and this is the only one that seemed to work
				(see index.php for related markup).
	*/
	$(progressBar).slider();
	$(progressBar).val(value);
	$(progressBar).attr('min', 0);
	$(progressBar).attr('max', 1);
	$(progressBar).attr('step', 1);
	$(progressBar).on('slidestart', function(event, ui){
		onSlideStart(canvas, progressBar);
	});
	$(progressBar).on('change', function(event, ui){
		onSlideChange(canvas, progressBar);
	});
	$(progressBar).on('slidestop', function(event, ui){
		onSlideStop(canvas, progressBar);
	});
	$(progressBar).slider('refresh'); //Slider must be "refreshed" after any JavaScript change to it, as it is an AJAX object.	

	playButton.click(function(){
		toggleMoviePlay(canvas);
	});
	
	canvas.progressBar = progressBar;
	canvas.playButton = playButton;
	canvas.timeLabel = timeText;
} //}}}
