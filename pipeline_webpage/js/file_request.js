var run_id = '';
var request_context = '';
var sample_id = '';
var chosen_report = '';

// async function launch_reporting() {

//     $("#report_selection").hide();
//     $("#report_selection_back").hide();
//     $("#email").hide();
//     $('#launch_report_button').hide();


//     var samples = document.getElementsByName('sample_id');
//     for(var i = 0; i < samples.length; i++) {
//         if(samples[i].checked==true){
//             sample_id.push(samples[i].id);
//         }
//     }

//     var email = document.getElementById('user_email').value;

//     var url = {};

//     var batch = new AWS.Batch({apiVersion: '2016-08-10'});
//     var s3 = new AWS.S3({apiVersion: '2006-03-01'});
//     var ses = new AWS.SES({apiVersion: '2010-12-01'});

//     for(var i = 0; i < sample_id.length; i++){

//         var panels = $("#"+sample_id[i]+"_select :selected").map((_, e) => e.value).get();

//         if(panels[0] == "low_coverage") {
//             var lc = "5x";
//             panels.splice(0, 1);
//         } else {
//             var lc = "30x";
//         }

//         url[sample_id[i]] = {};

//         var multiqc_url_params = { 
//             Bucket: 'hakmonkey-genetics-lab',
//             Key: 'Pipeline_Output/'+run_id+'/'+sample_id[i]+'/MultiQC/'+sample_id[i]+'.html',
//             Expires: 86400 // change to 86400 = 1 day
//         };

//         var multiqc_link = s3.getSignedUrl('getObject', multiqc_url_params);

//         url[sample_id[i]]['MultiQC'] = ['<a href='+multiqc_link+'>MultiQC Report</a>']

//         if(panels.length != 0){
//             for(var j = 0; j < panels.length; j++){
//                 var job_params = {
//                     jobDefinition: "var_class-ufl-germline:1", 
//                     jobName: sample_id[i]+'_'+panels[j], 
//                     jobQueue: "hakmonkey-var_class",
//                     containerOverrides: {
//                         'command': [
//                             'bash',
//                             '-c',
//                             'aws s3 cp s3://hakmonkey-genetics-lab/Pipeline_Output/'+run_id+'/'+sample_id[i]+'/variants/'+sample_id[i]+'_concat.vcf.gz /; aws s3 cp s3://hakmonkey-genetics-lab/Pipeline_Output/'+run_id+'/'+sample_id[i]+'/variants/'+sample_id[i]+'_concat.vcf.gz.tbi /; aws s3 cp s3://hakmonkey-genetics-lab/Pipeline/Reference/panels/'+panels[j]+' /; /varClass.py -v '+sample_id[i]+'_concat.vcf.gz -t 8 -s '+sample_id[i]+' -p '+panels[j]+' -c '+lc+'; /json_to_csv.py -j '+sample_id[i]+'_'+panels[j]+'_report.json; /g_ranges.py -j '+sample_id[i]+'_'+panels[j]+'_report.json -s '+sample_id[i]+'; /CNV_json_plot.R '+sample_id[i]+'; aws s3 cp '+sample_id[i]+'_'+panels[j]+'_report.json s3://hakmonkey-genetics-lab/Pipeline_Output/'+run_id+'/'+sample_id[i]+'/'+panels[j]+'/; aws s3 cp '+sample_id[i]+'_'+panels[j]+'_report.xlsx s3://hakmonkey-genetics-lab/Pipeline_Output/'+run_id+'/'+sample_id[i]+'/'+panels[j]+'/; aws s3 cp '+sample_id[i]+'_cnv.pdf s3://hakmonkey-genetics-lab/Pipeline_Output/'+run_id+'/'+sample_id[i]+'/'+panels[j]+'/'
//                         ]
//                     }
//                 };

//                 var json_url_params = { 
//                     Bucket: 'hakmonkey-genetics-lab',
//                     Key: 'Pipeline_Output/'+run_id+'/'+sample_id[i]+'/'+panels[j]+'/'+sample_id[i]+'_'+panels[j]+'_report.json',
//                     Expires: 86400 // change to 86400 = 1 day
//                 };
//                 var xlsx_url_params = { 
//                     Bucket: 'hakmonkey-genetics-lab',
//                     Key: 'Pipeline_Output/'+run_id+'/'+sample_id[i]+'/'+panels[j]+'/'+sample_id[i]+'_'+panels[j]+'_report.xlsx',
//                     Expires: 86400 // change to 86400 = 1 day
//                 };
//                 var cnv_url_params = { 
//                     Bucket: 'hakmonkey-genetics-lab',
//                     Key: 'Pipeline_Output/'+run_id+'/'+sample_id[i]+'/'+panels[j]+'/'+sample_id[i]+'_cnv.pdf',
//                     Expires: 86400 // change to 86400 = 1 day
//                 };
//                 var json_link = s3.getSignedUrl('getObject', json_url_params);
//                 var xlsx_link = s3.getSignedUrl('getObject', xlsx_url_params);
//                 var cnv_link = s3.getSignedUrl('getObject', cnv_url_params);
                
//                 url[sample_id[i]][panels[j]] = [
//                     '<a href='+json_link+'>JSON Report</a>',
//                     '<a href='+xlsx_link+'>XLSX Report</a>',
//                     '<a href='+cnv_link+'>CNV Plot</a>'
//                 ];

//                 // batch.submitJob(job_params, function(err, data) {
//                 //     if(err) {
//                 //         console.log(err, err.stack);
//                 //     } else {
//                 //         $("#launch_img").show();
//                 //         console.log(data);
//                 //     }
//                 // });
//             }
//         } else {
//             var job_params = {
//                 jobDefinition: "var_class-ufl-germline:1", 
//                 jobName: sample_id[i]+'_General_Report', 
//                 jobQueue: "hakmonkey-var_class",
//                 containerOverrides: {
//                     'command': [
//                         'bash',
//                         '-c',
//                         'aws s3 cp s3://hakmonkey-genetics-lab/Pipeline_Output/'+run_id+'/'+sample_id[i]+'/variants/'+sample_id[i]+'_concat.vcf.gz /; aws s3 cp s3://hakmonkey-genetics-lab/Pipeline_Output/'+run_id+'/'+sample_id[i]+'/variants/'+sample_id[i]+'_concat.vcf.gz.tbi /; /varClass.py -v '+sample_id[i]+'_concat.vcf.gz -t 8 -s '+sample_id[i]+' -c '+lc+'; /json_to_csv.py -j '+sample_id[i]+'_report.json; /g_ranges.py -j '+sample_id[i]+'_report.json -s '+sample_id[i]+'; /CNV_json_plot.R '+sample_id[i]+'; aws s3 cp '+sample_id[i]+'_report.json s3://hakmonkey-genetics-lab/Pipeline_Output/'+run_id+'/'+sample_id[i]+'/General_Report/; aws s3 cp '+sample_id[i]+'_report.xlsx s3://hakmonkey-genetics-lab/Pipeline_Output/'+run_id+'/'+sample_id[i]+'/General_Report/; aws s3 cp '+sample_id[i]+'_cnv.pdf s3://hakmonkey-genetics-lab/Pipeline_Output/'+run_id+'/'+sample_id[i]+'/General_Report/'
//                     ]
//                 }
//             };
            
//             var json_url_params = { 
//                 Bucket: 'hakmonkey-genetics-lab',
//                 Key: 'Pipeline_Output/'+run_id+'/'+sample_id[i]+'/General_Report/'+sample_id[i]+'_report.json',
//                 Expires: 86400 // change to 86400 = 1 day
//             };
//             var xlsx_url_params = { 
//                 Bucket: 'hakmonkey-genetics-lab',
//                 Key: 'Pipeline_Output/'+run_id+'/'+sample_id[i]+'/General_Report/'+sample_id[i]+'_report.xlsx',
//                 Expires: 86400 // change to 86400 = 1 day
//             };
//             var cnv_url_params = { 
//                 Bucket: 'hakmonkey-genetics-lab',
//                 Key: 'Pipeline_Output/'+run_id+'/'+sample_id[i]+'/General_Report/'+sample_id[i]+'_cnv.pdf',
//                 Expires: 86400 // change to 86400 = 1 day
//             };
//             var json_link = s3.getSignedUrl('getObject', json_url_params);
//             var xlsx_link = s3.getSignedUrl('getObject', xlsx_url_params);
//             var cnv_link = s3.getSignedUrl('getObject', cnv_url_params);

//             url[sample_id[i]]['General_Report'] = [
//                 '<a href='+json_link+'>JSON Report</a>',
//                 '<a href='+xlsx_link+'>XLSX Report</a>',
//                 '<a href='+cnv_link+'>CNV Plot</a>'
//             ];

//             // batch.submitJob(job_params, function(err, data) {
//             //     if(err) {
//             //         console.log(err, err.stack);
//             //     } else {
//             //         $("#launch_img").show();
//             //         console.log(data);
//             //     }
//             // });
//         }
//     }

//     console.log(email)

//     var url_list = JSON.stringify(url, null, 4);

//     console.log(url_list);

//     var email_params = {
//         Destination: {
//             ToAddresses: [ email ]
//         },
//         Message: {
//             Body: {
//                 Html: {
//                     Charset: "UTF-8", 
//                     Data: "The following is/are the link(s) for the requested report(s). The links should be valid for one day.<br/><br/>Please wait roughly 10 min before checking for the reports.<br/><br/><pre>"+url_list+"</pre><br/><br/>Cheers,<br/>Johnny"
//                 }
//             },
//             Subject: {
//                 Charset: "UTF-8",
//                 Data: "Link(s) for requested research report(s)"
//             }
//         },
//         Source: "jonathan.bravo@neurology.ufl.edu"
//     };

//     ses.sendEmail(email_params, function(err, data) {
//         if (err) console.log(err, err.stack); // an error occurred
//         else     console.log(data);           // successful response
//     });

//     $("#loader").show();

//     await new Promise(r => setTimeout(r, 2000));

//     $("#loader").hide();

//     $("#menu").show();
//     location.reload();
// }

// function get_panels(div_id, in_id){

//     var parent = document.getElementById(div_id);
//     var sample_box = document.getElementById(in_id);

//     if(sample_box.checked == true) {

//         var select_div = document.createElement("div");
//         select_div.setAttribute("class", "flex-child-child");
//         select_div.setAttribute("id", in_id+"_select-div");
//         var s = document.createElement("script");
//         var panels_select = document.createElement("select");

//         panels_select.setAttribute("multiple", "multiple");
//         panels_select.setAttribute("id", in_id+"_select")
//         s.setAttribute("type", "text/javascript");
//         s.setAttribute("id", in_id+"_my-script");
//         s.innerHTML = "$('#"+in_id+"_select').multiSelect();"

//         var low_coverage = document.createElement("option");
//         low_coverage.setAttribute("value", "low_coverage");
//         low_coverage.setAttribute("id", "low_coverage");
//         low_coverage.innerHTML = "LOW COVERAGE";
//         panels_select.appendChild(low_coverage);

//         for (const i in filtered_panels) {
//             var choiceSelection = document.createElement("option");

//             choiceSelection.setAttribute("value", filtered_panels[i]);
//             choiceSelection.setAttribute("id", "panel_"+i);

//             choiceSelection.innerHTML=filtered_panels[i];

//             panels_select.appendChild(choiceSelection);
//         }

//         select_div.appendChild(panels_select);
//         select_div.appendChild(s);

//         parent.appendChild(select_div);
//     } else {

//         var select_div = document.getElementById(in_id+"_select-div");

//         parent.removeChild(select_div);
//     }
// }

async function retrieve_report() {

    var request_email = document.getElementById("request_user_email").value;

    var ses = new AWS.SES({apiVersion: '2010-12-01'});

    // if(chosen_report == "multiqc_run") {

    // } else {

    // }

    // var cnv_link = s3.getSignedUrl('getObject', cnv_url_params);
    
    //     '<a href='+_link+'>JSON Report</a>',


    // var json_url_params = { 
    //     Bucket: 'hakmonkey-genetics-lab',
    //     Key: 'Pipeline_Output/'+run_id+'/'+sample_id[i]+'/'+panels[j]+'/'+sample_id[i]+'_'+panels[j]+'_report.json',
    //     Expires: 86400 // change to 86400 = 1 day
    // };
}

async function get_request_email() {

    $("#request_reports_box").hide();

    var chosen_reports = document.getElementsByName('report_name');
    if(chosen_reports.length != 0) {
        for(var i = 0; i < chosen_reports.length; i++) {
            if(chosen_reports[i].checked==true){
                chosen_report = chosen_reports[i].id;
            }
        }
    } else {
        chosen_report = "multiqc_run";
    }

    console.log(chosen_report);

    $("#request_email_box").show();
}

async function get_chosen_request_sample() {

    $("#request_samples_box").hide();

    var chosen_sample = document.getElementsByName('sample_id');
    for(var i = 0; i < chosen_sample.length; i++) {
        if(chosen_sample[i].checked==true){
            sample_id = chosen_sample[i].id;
        }
    }

    s3 = new AWS.S3({apiVersion: '2006-03-01'});

    var bucketParams = {
        Bucket : 'hakmonkey-genetics-lab',
        Prefix: 'Pipeline_Output/'+run_id+'/'+sample_id+'/',
        Delimiter: '/'
    };

    var reports = [];
    var available_reports = [];

    s3.listObjects(bucketParams, function(err, data) {
        if (err) {
            console.log("Error", err);
        } else {
            for(const i in data['CommonPrefixes']) {
                reports.push(data['CommonPrefixes'][i]['Prefix'].substring(25,).split("/")[1]);
            }
            for(const i in reports) {
                console.log(reports[i]);
                if(
                    reports[i].startsWith('_') == false &&
                    reports[i] != "alignment" &&
                    reports[i] != "ExpansionHunter" &&
                    reports[i].startsWith("fastqc") == false &&
                    reports[i] != "Panelcn_MOPS" &&
                    reports[i] != "Strelka2" &&
                    reports[i] != "Trimmomatic" &&
                    reports[i] != "wgs_metrics" &&
                    reports[i] != "wes_metrics" 
                ) {
                    available_reports.push(reports[i].slice(0,-1));
                }
            }
        }
    });

    $("#loader").show();

    await new Promise(r => setTimeout(r, 1000));

    var reports_list = document.getElementById('request_reports');

    for (var i = 0; i < available_reports.length; i++) {
        var choiceSelection = document.createElement('input');
        var choiceLabel = document.createElement('label');

        choiceSelection.setAttribute("type", "radio");
        choiceSelection.setAttribute("name", "report_name");
        choiceSelection.setAttribute("id", available_reports[i]);

        choiceLabel.appendChild(choiceSelection);
        choiceLabel.innerHTML += '\t'+available_reports[i]+"<br/><br/>";

        reports_list.appendChild(choiceLabel);
    }


    console.log(available_reports);

    $("#loader").hide();
    $("#request_reports_box").show();

}

async function get_request_sample() {

    s3 = new AWS.S3({apiVersion: '2006-03-01'});

    var bucketParams = {
        Bucket : 'hakmonkey-genetics-lab',
        Prefix: 'Pipeline_Output/'+run_id+'/',
        Delimiter: '/'
    };
    
    var samples = [];
    var filtered_samples = [];

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

    $("#loader").show();

    await new Promise(r => setTimeout(r, 3000));

    var samples_list = document.getElementById('request_samples');

    for (var i = 0; i < filtered_samples.length; i++) {
        var choiceSelection = document.createElement('input');
        var choiceLabel = document.createElement('label');

        choiceSelection.setAttribute("type", "radio");
        choiceSelection.setAttribute("name", "sample_id");
        choiceSelection.setAttribute("id", filtered_samples[i]);

        choiceLabel.appendChild(choiceSelection);
        choiceLabel.innerHTML += "\t"+filtered_samples[i]+"<br/><br/>";

        samples_list.appendChild(choiceLabel);
    }

    $("#loader").hide();
    $('#request_samples_box').show();
}

function r_or_s(){

    $("#request_run_or_sample_box").hide();

    var request_contexts = document.getElementsByName('r_or_s');
    for(var i = 0; i < request_contexts.length; i++) {
        if(request_contexts[i].checked==true){
            request_context = request_contexts[i].id;
        }
    }

    console.log(request_context);

    if(request_context == "sample") {
        get_request_sample();
    } else {
        get_request_email();
    }
}

function switch_to_contexts() {

    $("#request_runs_box").hide();

    var chosen_run = document.getElementsByName('run_id');
    for(var i = 0; i < chosen_run.length; i++) {
        if(chosen_run[i].checked==true){
            run_id = chosen_run[i].id;
        }
    }

    console.log(run_id);

    $("#request_run_or_sample_box").show();
}

async function get_request_runs() {

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

    var runs_list = document.getElementById("request_runs_list");

    for (const i in filtered_samples) {
        var choiceSelection = document.createElement("input");
        var choiceLabel = document.createElement("label");

        choiceSelection.setAttribute("type", "radio");
        choiceSelection.setAttribute("name", "run_id");
        choiceSelection.setAttribute("id", filtered_samples[i]);

        choiceLabel.appendChild(choiceSelection);
        choiceLabel.innerHTML += "\t"+filtered_samples[i]+"<br/><br/>";

        runs_list.appendChild(choiceLabel);
    }

    $("#loader").hide();
    $("#request_runs_box").show();
}