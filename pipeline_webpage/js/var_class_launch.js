var run_id = '';
var sample_id = [];
var filtered_panels = [];
var filtered_samples = [];

function distinct(value, index, self) {
    return self.indexOf(value) === index;
}

async function launch_reporting() {

    $("#report_selection").hide();
    $("#report_selection_back").hide();
    $("#email").hide();
    $('#launch_report_button').hide();


    var samples = document.getElementsByName('sample_id');
    for(var i = 0; i < samples.length; i++) {
        if(samples[i].checked==true){
            sample_id.push(samples[i].id);
        }
    }

    var email = document.getElementById('user_email').value;

    var url = {};

    var batch = new AWS.Batch({apiVersion: '2016-08-10'});
    var s3 = new AWS.S3({apiVersion: '2006-03-01'});
    var ses = new AWS.SES({apiVersion: '2010-12-01'});

    for(var i = 0; i < sample_id.length; i++){

        var panels = $("#"+sample_id[i]+"_select :selected").map((_, e) => e.value).get();

        if(panels[0] == "low_coverage") {
            var lc = "5x";
            panels.splice(0, 1);
        } else {
            var lc = "30x";
        }

        url[sample_id[i]] = {};

        var multiqc_url_params = { 
            Bucket: 'hakmonkey-genetics-lab',
            Key: 'Pipeline_Output/'+run_id+'/'+sample_id[i]+'/MultiQC/'+sample_id[i]+'.html',
            Expires: 86400 // change to 86400 = 1 day
        };

        var multiqc_link = s3.getSignedUrl('getObject', multiqc_url_params);

        url[sample_id[i]]['MultiQC'] = ['<a href='+multiqc_link+'>MultiQC Report</a>']

        if(panels.length != 0){
            for(var j = 0; j < panels.length; j++){
                var job_params = {
                    jobDefinition: "var_class-ufl-germline:1", 
                    jobName: sample_id[i]+'_'+panels[j], 
                    jobQueue: "hakmonkey-var_class",
                    containerOverrides: {
                        'command': [
                            'bash',
                            '-c',
                            'aws s3 cp s3://hakmonkey-genetics-lab/Pipeline_Output/'+run_id+'/'+sample_id[i]+'/variants/'+sample_id[i]+'_concat.vcf.gz /; aws s3 cp s3://hakmonkey-genetics-lab/Pipeline_Output/'+run_id+'/'+sample_id[i]+'/variants/'+sample_id[i]+'_concat.vcf.gz.tbi /; aws s3 cp s3://hakmonkey-genetics-lab/Pipeline/Reference/panels/'+panels[j]+' /; aws s3 sync s3://hakmonkey-genetics-lab/Pipeline/Reporting/ /; /reporting.py -v '+sample_id[i]+'_concat.vcf.gz -t 16 -s '+sample_id[i]+' -p '+panels[j]+' -c '+lc+'; /json_to_csv.py -j '+sample_id[i]+'_'+panels[j]+'_report.json; /g_ranges.py -j '+sample_id[i]+'_'+panels[j]+'_report.json -s '+sample_id[i]+'; /CNV_json_plot.R '+sample_id[i]+'; aws s3 cp '+sample_id[i]+'_'+panels[j]+'_report.json s3://hakmonkey-genetics-lab/Pipeline_Output/'+run_id+'/'+sample_id[i]+'/'+panels[j]+'/; aws s3 cp '+sample_id[i]+'_'+panels[j]+'_report.xlsx s3://hakmonkey-genetics-lab/Pipeline_Output/'+run_id+'/'+sample_id[i]+'/'+panels[j]+'/; aws s3 cp '+sample_id[i]+'_cnv.pdf s3://hakmonkey-genetics-lab/Pipeline_Output/'+run_id+'/'+sample_id[i]+'/'+panels[j]+'/'
                        ]
                    }
                };

                var json_url_params = { 
                    Bucket: 'hakmonkey-genetics-lab',
                    Key: 'Pipeline_Output/'+run_id+'/'+sample_id[i]+'/'+panels[j]+'/'+sample_id[i]+'_'+panels[j]+'_report.json',
                    Expires: 86400 // change to 86400 = 1 day
                };
                var xlsx_url_params = { 
                    Bucket: 'hakmonkey-genetics-lab',
                    Key: 'Pipeline_Output/'+run_id+'/'+sample_id[i]+'/'+panels[j]+'/'+sample_id[i]+'_'+panels[j]+'_report.xlsx',
                    Expires: 86400 // change to 86400 = 1 day
                };
                var cnv_url_params = { 
                    Bucket: 'hakmonkey-genetics-lab',
                    Key: 'Pipeline_Output/'+run_id+'/'+sample_id[i]+'/'+panels[j]+'/'+sample_id[i]+'_cnv.pdf',
                    Expires: 86400 // change to 86400 = 1 day
                };
                var json_link = s3.getSignedUrl('getObject', json_url_params);
                var xlsx_link = s3.getSignedUrl('getObject', xlsx_url_params);
                var cnv_link = s3.getSignedUrl('getObject', cnv_url_params);
                
                url[sample_id[i]][panels[j]] = [
                    '<a href='+json_link+'>JSON Report</a>',
                    '<a href='+xlsx_link+'>XLSX Report</a>',
                    '<a href='+cnv_link+'>CNV Plot</a>'
                ];

                // batch.submitJob(job_params, function(err, data) {
                //     if(err) {
                //         console.log(err, err.stack);
                //     } else {
                //         $("#launch_img").show();
                //         console.log(data);
                //     }
                // });
            }
        } else {
            var job_params = {
                jobDefinition: "var_class-ufl-germline:1", 
                jobName: sample_id[i]+'_General_Report', 
                jobQueue: "hakmonkey-var_class",
                containerOverrides: {
                    'command': [
                        'bash',
                        '-c',
                        'aws s3 cp s3://hakmonkey-genetics-lab/Pipeline_Output/'+run_id+'/'+sample_id[i]+'/variants/'+sample_id[i]+'_concat.vcf.gz /; aws s3 cp s3://hakmonkey-genetics-lab/Pipeline_Output/'+run_id+'/'+sample_id[i]+'/variants/'+sample_id[i]+'_concat.vcf.gz.tbi /; aws s3 sync s3://hakmonkey-genetics-lab/Pipeline/Reporting/ /; /reporting.py -v '+sample_id[i]+'_concat.vcf.gz -t 8 -s '+sample_id[i]+' -c '+lc+'; /json_to_csv.py -j '+sample_id[i]+'_report.json; /g_ranges.py -j '+sample_id[i]+'_report.json -s '+sample_id[i]+'; /CNV_json_plot.R '+sample_id[i]+'; aws s3 cp '+sample_id[i]+'_report.json s3://hakmonkey-genetics-lab/Pipeline_Output/'+run_id+'/'+sample_id[i]+'/General_Report/; aws s3 cp '+sample_id[i]+'_report.xlsx s3://hakmonkey-genetics-lab/Pipeline_Output/'+run_id+'/'+sample_id[i]+'/General_Report/; aws s3 cp '+sample_id[i]+'_cnv.pdf s3://hakmonkey-genetics-lab/Pipeline_Output/'+run_id+'/'+sample_id[i]+'/General_Report/'
                    ]
                }
            };
            
            var json_url_params = { 
                Bucket: 'hakmonkey-genetics-lab',
                Key: 'Pipeline_Output/'+run_id+'/'+sample_id[i]+'/General_Report/'+sample_id[i]+'_report.json',
                Expires: 86400 // change to 86400 = 1 day
            };
            var xlsx_url_params = { 
                Bucket: 'hakmonkey-genetics-lab',
                Key: 'Pipeline_Output/'+run_id+'/'+sample_id[i]+'/General_Report/'+sample_id[i]+'_report.xlsx',
                Expires: 86400 // change to 86400 = 1 day
            };
            var cnv_url_params = { 
                Bucket: 'hakmonkey-genetics-lab',
                Key: 'Pipeline_Output/'+run_id+'/'+sample_id[i]+'/General_Report/'+sample_id[i]+'_cnv.pdf',
                Expires: 86400 // change to 86400 = 1 day
            };
            var json_link = s3.getSignedUrl('getObject', json_url_params);
            var xlsx_link = s3.getSignedUrl('getObject', xlsx_url_params);
            var cnv_link = s3.getSignedUrl('getObject', cnv_url_params);

            url[sample_id[i]]['General_Report'] = [
                '<a href='+json_link+'>JSON Report</a>',
                '<a href='+xlsx_link+'>XLSX Report</a>',
                '<a href='+cnv_link+'>CNV Plot</a>'
            ];

            // batch.submitJob(job_params, function(err, data) {
            //     if(err) {
            //         console.log(err, err.stack);
            //     } else {
            //         $("#launch_img").show();
            //         console.log(data);
            //     }
            // });
        }
    }

    console.log(email)

    var url_list = JSON.stringify(url, null, 4);

    console.log(url_list);

    var email_params = {
        Destination: {
            ToAddresses: [ email ]
        },
        Message: {
            Body: {
                Html: {
                    Charset: "UTF-8", 
                    Data: "The following is/are the link(s) for the requested report(s). The links should be valid for one day.<br/><br/>Please wait roughly 10 min before checking for the reports.<br/><br/><pre>"+url_list+"</pre><br/><br/>Cheers,<br/>Johnny"
                }
            },
            Subject: {
                Charset: "UTF-8",
                Data: "Link(s) for requested research report(s)"
            }
        },
        Source: "jonathan.bravo@neurology.ufl.edu"
    };

    ses.sendEmail(email_params, function(err, data) {
        if (err) console.log(err, err.stack); // an error occurred
        else     console.log(data);           // successful response
    });

    $("#loader").show();

    await new Promise(r => setTimeout(r, 2000));

    $("#loader").hide();

    $("#menu").show();
    location.reload();
}

function select_all_samples() {

    for(var k = 0; k < filtered_samples.length; k++){
        var sample_check = document.getElementById(String(filtered_samples[k]));
        // sample_check.checked = !sample_check.checked;
        if (sample_check.checked == false) {
            sample_check.checked = true;

            if ("createEvent" in document) {
                var evt = document.createEvent("HTMLEvents");
                evt.initEvent("change", false, true);
                sample_check.dispatchEvent(evt);
            } else {
                sample_check.fireEvent("onchange");
            }
        }
    }
}

function get_panels_parent(){

    var ms_filtered_samples = [];

    for(var i = 0; i < filtered_samples.length; i++){
        ms_filtered_samples.push("'ms-"+filtered_samples[i]+"_select'");
    }

    var select_all_box = document.createElement("div");
    select_all_box.setAttribute("class", "flex-child-child");
    select_all_box.setAttribute("id", "select_all_box");

    var select_all = document.createElement("button");
    select_all.setAttribute("onclick", "select_all_samples()");
    select_all.setAttribute("id", "select_all_button");
    select_all.innerHTML = "Select All";

    select_all_box.appendChild(select_all);

    var parent_box = document.getElementById("parent_selector_box");
    var child_box = document.getElementById("parent_selector_child_box");
    var parent_selector = document.getElementById("parent_selector");

    var low_coverage = document.createElement("option");
    low_coverage.setAttribute("value", "5x WGS");
    low_coverage.setAttribute("id", "low_coverage");
    low_coverage.innerHTML = "5x WGS";
    parent_selector.appendChild(low_coverage);

    for (const i in filtered_panels) {
        var choiceSelection = document.createElement("option");

        choiceSelection.setAttribute("value", filtered_panels[i]);
        choiceSelection.setAttribute("id", "panel_"+i);

        choiceSelection.innerHTML=filtered_panels[i];

        parent_selector.appendChild(choiceSelection);
    }

    var select_code = document.createElement("script");
    select_code.setAttribute("type", "text/javascript");
    select_code.setAttribute("id", "parent_selector_script");
    select_code.innerHTML = `
        $('#parent_selector').multiSelect({
            selectableHeader: "<input type='text' class='search-input' autocomplete='off' placeholder='panels to select'>",
            selectionHeader: "<input type='text' class='search-input' autocomplete='off' placeholder='selected panels'>",
            afterInit: function(ms){
                var that = this,
                    $selectableSearch = that.$selectableUl.prev(),
                    $selectionSearch = that.$selectionUl.prev(),
                    selectableSearchString = '#'+that.$container.attr('id')+' .ms-elem-selectable:not(.ms-selected)',
                    selectionSearchString = '#'+that.$container.attr('id')+' .ms-elem-selection.ms-selected';

                that.qs1 = $selectableSearch.quicksearch(selectableSearchString)
                .on('keydown', function(e){
                    if (e.which === 40){
                        that.$selectableUl.focus();
                        return false;
                    }
                });

                that.qs2 = $selectionSearch.quicksearch(selectionSearchString)
                .on('keydown', function(e){
                    if (e.which == 40){
                        that.$selectionUl.focus();
                        return false;
                    }
                });
            },
            afterSelect: function(ms){
                this.qs1.cache();
                this.qs2.cache();

                var selectors = [${ms_filtered_samples}];

                for(var j = 0; j < selectors.length; j++){
                    var p = document.getElementById(String(selectors[j])).children; //selector
                    var s = p[0].children[0]; // selectable list
                    var sl = s.children; // listed elements
                    var l = p[1].children[0];
                    var ll = l.children;

                    for(var i = 0; i< sl.length; i++) {
                        var v = sl[i].children[0].innerHTML; // value
                        if(v == ms[0]){
                            sl[i].style = "display: none;";
                            ll[i].className = "ms-elem-selection ms-selected";
                            ll[i].style = "";
                        }
                    }
                }
            },
            afterDeselect: function(ms){
                this.qs1.cache();
                this.qs2.cache();

                var selectors = [${ms_filtered_samples}];

                for(var j = 0; j < selectors.length; j++){
                    var p = document.getElementById(String(selectors[j])).children; //selector
                    var s = p[0].children[0]; // selectable list
                    var sl = s.children; // listed elements
                    var l = p[1].children[0];
                    var ll = l.children;

                    for(var i = 0; i< ll.length; i++) {
                        var v = ll[i].children[0].innerHTML; // value
                        if(v == ms[0]){
                            ll[i].style = "display: none;";
                            ll[i].className = "ms-elem-selection";
                            sl[i].style = "";
                        }
                    }
                }
            }
        })`;
        child_box.appendChild(select_code);
        parent_box.appendChild(select_all_box);
}

function get_panels(div_id, in_id){

    var parent = document.getElementById(div_id);
    var sample_box = document.getElementById(in_id);

    if(sample_box.checked == true) {

        var select_div = document.createElement("div");
        select_div.setAttribute("class", "flex-child-child");
        select_div.setAttribute("id", in_id+"_select-div");
        var s = document.createElement("script");
        var panels_select = document.createElement("select");

        panels_select.setAttribute("multiple", "multiple");
        panels_select.setAttribute("id", in_id+"_select")
        s.setAttribute("type", "text/javascript");
        s.setAttribute("id", in_id+"_my-script");
        s.innerHTML = "$('#"+in_id+"_select').multiSelect();"

        var low_coverage = document.createElement("option");
        low_coverage.setAttribute("value", "low_coverage");
        low_coverage.setAttribute("id", "low_coverage");
        low_coverage.innerHTML = "5x WGS";
        panels_select.appendChild(low_coverage);

        for (const i in filtered_panels) {
            var choiceSelection = document.createElement("option");

            choiceSelection.setAttribute("value", filtered_panels[i]);
            choiceSelection.setAttribute("id", "panel_"+i);

            choiceSelection.innerHTML=filtered_panels[i];

            panels_select.appendChild(choiceSelection);
        }

        select_div.appendChild(panels_select);
        select_div.appendChild(s);

        parent.appendChild(select_div);
    } else {

        var select_div = document.getElementById(in_id+"_select-div");

        parent.removeChild(select_div);
    }
}

async function get_report_sample() {

    $('#report_runs_box').hide();
    

    var runs = document.getElementsByName('run_id');
    for(var i in runs) {
        if(runs[i].checked==true){
            run_id = runs[i].id;
        }
    }

    s3 = new AWS.S3({apiVersion: '2006-03-01'});

    //SAMPLES

    var bucketParams = {
        Bucket : 'hakmonkey-genetics-lab',
        Prefix: 'Pipeline_Output/'+run_id+'/',
        Delimiter: '/'
    };
    
    var samples = [];
    filtered_samples = [];

    s3.listObjects(bucketParams, function(err, data) {
        if (err) {
            console.log("Error", err);
        } else {
            for(const i in data['CommonPrefixes']) {
                samples.push(data['CommonPrefixes'][i]['Prefix'].substring(25,));
            }
            samples = samples.filter(distinct);
            samples = samples.filter(item => item);
            for(const i in samples) {
                if(samples[i].startsWith('_') == false && samples[i] != 'MultiQC/') {
                    filtered_samples.push(samples[i].slice(0,-1));
                }
            }
        }
    });

    // PANELS

    var bucketParams = {
        Bucket : 'hakmonkey-genetics-lab',
        Prefix: 'Pipeline/Reference/panels/',
        Delimiter: '/'
    };
    
    var panels = [];
    filtered_panels = [];

    s3.listObjects(bucketParams, function(err, data) {
        if (err) {
            console.log("Error", err);
        } else {
            for(const i in data['Contents']) {
                panels.push(data['Contents'][i]['Key'].substring(26,));
            }
            panels = panels.filter(distinct);
            panels = panels.filter(item => item);
            for(const i in panels) {
                if(panels[i].startsWith('_') == false) {
                    filtered_panels.push(panels[i]);
                }
            }
        }
    });

    $("#loader").show();

    await new Promise(r => setTimeout(r, 3000));

    var sample_list = document.getElementById('samples');

    for (const i in filtered_samples) {
        var sample_div = document.createElement('div');
        sample_div.setAttribute("class", "flex-child");
        sample_div.setAttribute("id", filtered_samples[i]+"_div");

        var sample_child_div = document.createElement("div");
        sample_child_div.setAttribute("class","flex-child-child");

        //var choiceSelection = document.createElement('input');
        //var choiceLabel = document.createElement('label');

        //choiceSelection.setAttribute("type", "checkbox");
        //choiceSelection.setAttribute("name", "sample_id");
        //choiceSelection.setAttribute("id", filtered_samples[i]);
        //choiceSelection.setAttribute("onchange", "get_panels('"+filtered_samples[i]+"_div','"+filtered_samples[i]+"')");

        //choiceLabel.appendChild(choiceSelection);
        //choiceLabel.innerHTML += '\t'+filtered_samples[i]+"<br/><br/>";

        // sample_child_div.appendChild(choiceLabel);
        // sample_div.appendChild(sample_child_div);
        // sample_list.appendChild(sample_div);

        var checkboxClass = document.createElement('label');
        checkboxClass.setAttribute('class', 'checkbox');

        var checkboxInput = document.createElement('span');
        checkboxInput.setAttribute('class', 'checkbox__input');

        var choiceSelection = document.createElement('input');
        choiceSelection.setAttribute("type", "checkbox");
        choiceSelection.setAttribute("name", "sample_id");
        choiceSelection.setAttribute("id", filtered_samples[i]);
        choiceSelection.setAttribute("onchange", "get_panels('"+filtered_samples[i]+"_div','"+filtered_samples[i]+"')");

        var checkboxControl = document.createElement('span');
        checkboxControl.setAttribute('class', 'checkbox__control');

        var checkboxSvg = document.createElement('svg');
        checkboxSvg.setAttribute('xmlns', 'http://www.w3.org/2000/svg');
        checkboxSvg.setAttribute('viewBox', '0 0 24 24');
        checkboxSvg.setAttribute('aria-hidden', 'true');
        checkboxSvg.setAttribute('focusable', 'false');

        var checkboxPath = document.createElement('path');
        checkboxPath.setAttribute('fill', 'none');
        checkboxPath.setAttribute('stroke', 'currentColor');
        checkboxPath.setAttribute('stroke-width', '3');
        checkboxPath.setAttribute('d', 'M1.73 12.91l6.37 6.37L22.79 4.59');

        var checkboxLabel = document.createElement('span');
        checkboxLabel.setAttribute('class', 'radio__label');
        checkboxLabel.innerHTML += '\t'+filtered_samples[i];

        checkboxSvg.appendChild(checkboxPath);
        checkboxControl.appendChild(checkboxSvg);

        checkboxInput.appendChild(choiceSelection);
        checkboxInput.appendChild(checkboxControl);

        checkboxClass.appendChild(checkboxInput);
        checkboxClass.appendChild(checkboxLabel);

        sample_child_div.appendChild(checkboxClass);
        sample_div.appendChild(sample_child_div);
        sample_list.appendChild(sample_div);
    }

    get_panels_parent();

    $("#loader").hide();
    $("#report_selection_back").show();
    $('#report_selection').show();
    $('#report_sample_box').show();
    $('#email').show();
    $('#launch_report_button').show();

}

async function get_report_runs() {

    s3 = new AWS.S3({apiVersion: '2006-03-01'});

    var bucketParams = {
        Bucket : 'hakmonkey-genetics-lab',
        Prefix: 'Pipeline_Output/',
        Delimiter: '/'
    };
    
    var samples = [];
    var filtered_samples = [];

    s3.listObjects(bucketParams, function(err, data) {
        if (err) {
            console.log("Error", err);
        } else {
            for(const i in data['CommonPrefixes']) {
                samples.push(data['CommonPrefixes'][i]['Prefix'].substring(16, 24));
            }
            samples = samples.filter(distinct);
            samples = samples.filter(item => item);
            for(const i in samples) {
                if(samples[i].startsWith('_') == false) {
                    filtered_samples.push(samples[i]);
                }
            }
        }
    });

    $("#loader").show();

    await new Promise(r => setTimeout(r, 2000));

    var runs_list = document.getElementById('report_runs_list');

    for (const i in filtered_samples) {
        var radioClass = document.createElement('label');
        var radioInput = document.createElement('span');
        var choiceSelection = document.createElement('input');
        var radioControl = document.createElement('span');
        var radioLabel = document.createElement('span');

        radioClass.setAttribute('class', 'radio');
        radioInput.setAttribute('class', 'radio__input');
        radioControl.setAttribute('class', 'radio__control');
        radioLabel.setAttribute('class', 'radio__label');

        choiceSelection.setAttribute('type', 'radio');
        choiceSelection.setAttribute('name', 'run_id');
        choiceSelection.setAttribute('id', filtered_samples[i]);

        radioLabel.innerHTML += '\t'+filtered_samples[i];

        radioInput.appendChild(choiceSelection);
        radioInput.appendChild(radioControl);

        radioClass.appendChild(radioInput);
        radioClass.appendChild(radioLabel);

        runs_list.appendChild(radioClass); 
    }

    $("#loader").hide();
    $("#report_runs_box").show();
}