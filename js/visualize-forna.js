
	// global FORNA container to be used for visualization
	fornaStructureContainer = null;

	
    /* *************************************************************************
     * This function initialize the forna container.
     * *************************************************************************/
    function initalLoadOfFornaContainer() {
    	if (typeof(fornaStructureContainer) === 'undefined' || fornaStructureContainer === null) {
    		fornaStructureContainer = new FornaContainer("#fornaInputStructure", {'applyForce': false});
        }
    }

    /* ********************************************************************
     * Clears the FORNA output
     */
    function fornaClear() {
    	// ensure container is set up
    	initalLoadOfFornaContainer();
    	// start visualization
    	// clear old data
    	fornaStructureContainer.clearNodes();
    }
    
    
    /* ********************************************************************
     * Renders a structure using FORNA
     */
    function fornaRendering( sequence, structureDomElement ) {
    	
    	var structure = structureDomElement.innerHTML;
    	
    	if (structure.indexOf('X')!=-1) {
    		structure = structure.replace(/X/g,'.');
    	}
    	
		// init visualization
		var options = {
			'structure' : structure,
			'sequence' : sequence,
			'name' : 'selected structure'
		};

    	// ensure container is set up
		fornaClear();
    	// add current structure
    	fornaStructureContainer.addRNA(options.structure, options);
    }
