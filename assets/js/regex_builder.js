window.onload = function() {
    
    link_checkboxes("fofchk", "fofchk_sub");
    link_checkboxes("rockchk", "rockchk_sub");
    
    // Register regex re-generator on checkbox click
    // Assumes all checkboxes on the page want to be registered
    var boxes = document.querySelectorAll('input[type=checkbox]')
    for(var i = 0; i < boxes.length; i++) {
        var box = boxes[i];
        box.addEventListener("change", build_regex);
    }
    
    build_regex();
}

function link_checkboxes(parent_id, child_name){
   // Register linking function for FoF checkboxes
    var subboxes = document.getElementsByName(child_name),
        parentbox = document.getElementById(parent_id);

    parentbox.addEventListener("change", function() {
                                             for(var i = 0; i < subboxes.length; i++)
                                                 subboxes[i].checked = this.checked;
                                         });
    for(var i=0; i<subboxes.length; i++) {
      subboxes[i].addEventListener("change", function() {
                                                 var checkedCount = document.querySelectorAll('[name="' + child_name + '"]:checked').length;
                                                 for(var i=0; i<subboxes.length; i++) {
                                                     parentbox.checked = checkedCount > 0;
                                                     parentbox.indeterminate = checkedCount > 0 && checkedCount < subboxes.length;
                                                 }
                                             });
    }
}

function build_regex() {
    //console.log('===============');
    // 'this' refers to the current checkbox in the EventListener model
    // but we'll probably be iterating over all checkboxes anyway
    var boxes = document.getElementsByClassName('chk');
    
    var base_command = "wget -r -nH -np --cut-dirs=2 -R 'index.html*' --accept-regex=";
    var base_url = "http://nbody.rc.fas.harvard.edu/public/AbacusCosmos/";
    var regex = "(_products/$|info/";  // always include _product and info directories
    var urls = "";
    var reject_regex = "";  // might be a no-op
    
    // Add URLs to the command based on which simulation sets
    var sims = document.getElementsByName("sims")[0];
    var sim_boxes = sims.getElementsByClassName("chk");
    for(var i = 0; i < sim_boxes.length; i++) {
        var box = sim_boxes[i];
        if(box.checked){
            var path = box.getAttribute("data-path");
            urls += base_url + path + "/ ";
        }
    }
    
    // Add products to the regex
    var products = document.getElementsByName("products")[0];
    var prod_boxes = products.getElementsByClassName("chk");
    for(var i = 0; i < prod_boxes.length; i++) {
        var box = prod_boxes[i];
        var product = box.getAttribute("data-product");
        //console.log('Product: ' + product + '; Checked: ' + box.checked);
        if(box.checked){
            //alert(product);
            regex += "|_" + product + "/$";
            
            // Add a reject regex if any product families have unselected files
            var sub_boxes = document.getElementsByName(box.id + "_sub");
            for(var j = 0; j < sub_boxes.length; j++) {
                var subbox = sub_boxes[j];
                if(!subbox.checked){
                    var fn = subbox.getAttribute("data-fn");
                    var pattern = product + '/.*/' + fn
                    // If this is the first reject regex, init the string
                    if(!reject_regex)
                        reject_regex += "--reject-regex='(" + pattern;
                    else
                        reject_regex += "|" + pattern;
                }
            }
        }
    }
    if(reject_regex)
        reject_regex += ")'"
    
    // Add redshifts to the regex
    var redshifts = document.getElementsByName("redshifts")[0];
    var z_boxes = redshifts.getElementsByClassName("chk");
    for(var i = 0; i < z_boxes.length; i++) {
        var box = z_boxes[i];
        if(box.checked){
            var redshift = box.getAttribute("data-redshift");
            regex += "|" + redshift + "/";
        }
    }
    
    regex += ")";
    
    // Construct the final command
    var command = base_command + "'" + regex + "' " + reject_regex + ' ' + urls;
    
    // and display it
    var textbox = document.getElementById("wget-command");
    textbox.value = command;
    //console.log('===============');
}
