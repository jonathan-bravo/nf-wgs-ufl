var run_id = '';
var sample_id = '';
var panel_choices = [];

function distinct(value, index, self) {
    return self.indexOf(value) === index;
}

function launch_reporting() {

    $("#report_selection").hide();
    $("#email").hide();
    $('#launch_report_button').hide();


    var samples = document.getElementsByName('sample_id');
    for(var i in samples) {
        if(samples[i].checked==true){
            sample_id = samples[i].id;
        }
    }

    var panels = document.getElementsByName('panel_name');
    for(var i in panels) {
        if(panels[i].checked==true){
            panel_choices.push(panels[i].id);
        }
    }

    var email = document.getElementById('user_email').value;

    var url = {};

    var batch = new AWS.Batch({apiVersion: '2016-08-10'});
    var s3 = new AWS.S3({apiVersion: '2006-03-01'});
    var ses = new AWS.SES({apiVersion: '2010-12-01'});


    if(panel_choices.length != 0){
        for(var i in panel_choices){
            var job_params = {
                jobDefinition: "var_class-ufl-germline:1", 
                jobName: sample_id+'_'+panel_choices[i], 
                jobQueue: "hakmonkey-var_class",
                containerOverrides: {
                    'command': [
                        'bash',
                        '-c',
                        'aws s3 cp s3://hakmonkey-genetics-lab/Pipeline_Output/'+run_id+'/'+sample_id+'/variants/'+sample_id+'_concat.vcf.gz /; aws s3 cp s3://hakmonkey-genetics-lab/Pipeline_Output/'+run_id+'/'+sample_id+'/variants/'+sample_id+'_concat.vcf.gz.tbi /; aws s3 cp s3://hakmonkey-genetics-lab/Pipeline/Reference/panels/'+panel_choices[i]+' /; /varClass.py -v '+sample_id+'_concat.vcf.gz -t 8 -s '+sample_id+' -p '+panel_choices[i]+'; /json_to_csv.py -j '+sample_id+'_'+panel_choices[i]+'_report.json; /g_ranges.py -j '+sample_id+'_'+panel_choices[i]+'_report.json -s '+sample_id+'; /CNV_json_plot.R '+sample_id+'; aws s3 cp '+sample_id+'_'+panel_choices[i]+'_report.json s3://hakmonkey-genetics-lab/Pipeline_Output/'+run_id+'/'+sample_id+'/'+panel_choices[i]+'/; aws s3 cp '+sample_id+'_'+panel_choices[i]+'_report.xlsx s3://hakmonkey-genetics-lab/Pipeline_Output/'+run_id+'/'+sample_id+'/'+panel_choices[i]+'/; aws s3 cp '+sample_id+'_cnv.pdf s3://hakmonkey-genetics-lab/Pipeline_Output/'+run_id+'/'+sample_id+'/'+panel_choices[i]+'/'
                    ]
                }
            };

            var json_url_params = { 
                Bucket: 'hakmonkey-genetics-lab',
                Key: 'Pipeline_Output/'+run_id+'/'+sample_id+'/'+panel_choices[i]+'/'+sample_id+'_'+panel_choices[i]+'_report.json',
                Expires: 86400 // change to 86400 = 1 day
            };
            var xlsx_url_params = { 
                Bucket: 'hakmonkey-genetics-lab',
                Key: 'Pipeline_Output/'+run_id+'/'+sample_id+'/'+panel_choices[i]+'/'+sample_id+'_'+panel_choices[i]+'_report.xlsx',
                Expires: 86400 // change to 86400 = 1 day
            };
            var cnv_url_params = { 
                Bucket: 'hakmonkey-genetics-lab',
                Key: 'Pipeline_Output/'+run_id+'/'+sample_id+'/'+panel_choices[i]+'/'+sample_id+'_cnv.pdf',
                Expires: 86400 // change to 86400 = 1 day
            };
            var json_link = s3.getSignedUrl('getObject', json_url_params);
            var xlsx_link = s3.getSignedUrl('getObject', xlsx_url_params);
            var cnv_link = s3.getSignedUrl('getObject', cnv_url_params);
            
            url[panel_choices[i]] = [
                '<a href='+json_link+'>JSON Report</a>',
                '<a href='+xlsx_link+'>XLSX Report</a>',
                '<a href='+cnv_link+'>CNV Plot</a>'
            ];

            batch.submitJob(job_params, function(err, data) {
                if(err) {
                    console.log(err, err.stack);
                } else {
                    $("#launch_img").show();
                    console.log(data);
                }
            });
            
        }
    } else {
        var job_params = {
            jobDefinition: "var_class-ufl-germline:1", 
            jobName: sample_id+'_'+panel_choices[i], 
            jobQueue: "hakmonkey-var_class",
            containerOverrides: {
                'command': [
                    'bash',
                    '-c',
                    'aws s3 cp s3://hakmonkey-genetics-lab/Pipeline_Output/'+run_id+'/'+sample_id+'/variants/'+sample_id+'_concat.vcf.gz /; aws s3 cp s3://hakmonkey-genetics-lab/Pipeline_Output/'+run_id+'/'+sample_id+'/variants/'+sample_id+'_concat.vcf.gz.tbi /; /varClass.py -v '+sample_id+'_concat.vcf.gz -t 8 -s '+sample_id+'; /json_to_csv.py -j '+sample_id+'_report.json; /g_ranges.py -j '+sample_id+'_report.json -s '+sample_id+'; /CNV_json_plot.R '+sample_id+'; aws s3 cp '+sample_id+'_report.json s3://hakmonkey-genetics-lab/Pipeline_Output/'+run_id+'/'+sample_id+'/General_Report/; aws s3 cp '+sample_id+'_report.xlsx s3://hakmonkey-genetics-lab/Pipeline_Output/'+run_id+'/'+sample_id+'/General_Report/; aws s3 cp '+sample_id+'_cnv.pdf s3://hakmonkey-genetics-lab/Pipeline_Output/'+run_id+'/'+sample_id+'/General_Report/'
                ]
            }
        };
        
        var json_url_params = { 
            Bucket: 'hakmonkey-genetics-lab',
            Key: 'Pipeline_Output/'+run_id+'/'+sample_id+'/General_Report/'+sample_id+'_report.json',
            Expires: 86400 // change to 86400 = 1 day
        };
        var xlsx_url_params = { 
            Bucket: 'hakmonkey-genetics-lab',
            Key: 'Pipeline_Output/'+run_id+'/'+sample_id+'/General_Report/'+sample_id+'_report.xlsx',
            Expires: 86400 // change to 86400 = 1 day
        };
        var cnv_url_params = { 
            Bucket: 'hakmonkey-genetics-lab',
            Key: 'Pipeline_Output/'+run_id+'/'+sample_id+'/General_Report/'+sample_id+'_cnv.pdf',
            Expires: 86400 // change to 86400 = 1 day
        };
        var json_link = s3.getSignedUrl('getObject', json_url_params);
        var xlsx_link = s3.getSignedUrl('getObject', xlsx_url_params);
        var cnv_link = s3.getSignedUrl('getObject', cnv_url_params);

        url['General_Report'] = [
            '<a href='+json_link+'>JSON Report</a>',
            '<a href='+xlsx_link+'>XLSX Report</a>',
            '<a href='+cnv_link+'>CNV Plot</a>'
        ];

        batch.submitJob(job_params, function(err, data) {
            if(err) {
                console.log(err, err.stack);
            } else {
                $("#launch_img").show();
                console.log(data);
            }
        });

    }

    var url_list = JSON.stringify(url, null, 4);

    var email_params = {
        Destination: {
            ToAddresses: [ email ]
        },
        Message: {
            Body: {
                Html: {
                    Charset: "UTF-8", 
                    Data: "The following is/are the link(s) for the requested report(s).The links should be valid for one day.<br/><br/>Please wait at least 10-15 min before checking for the reports.<br/><br/><pre>"+url_list+"</pre><br/><br/>Cheers,<br/>Johnny"
                }
            },
            Subject: {
                Charset: "UTF-8",
                Data: "Link for "+sample_id+" research report"
            }
        },
        Source: "jonathan.bravo@neurology.ufl.edu"
    };

    ses.sendEmail(email_params, function(err, data) {
        if (err) console.log(err, err.stack); // an error occurred
        else     console.log(data);           // successful response
    });
    
    
}


async function get_report_sample() {

    $('#report_runs_box').hide();
    $('#report_selection').show();
    

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

    // PANELS

    var bucketParams = {
        Bucket : 'hakmonkey-genetics-lab',
        Prefix: 'Pipeline/Reference/panels/',
        Delimiter: '/'
    };
    
    var panels = [];
    var filtered_panels = [];

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

    $('#report_sample_box').show();
    $('#report_panel_box').show();

    var samples_box = document.getElementById("report_sample_box");
    var sample_list = document.getElementById('samples');

    var panels_box = document.getElementById('report_panel_box');
    var panel_list = document.getElementById('panels');

    for (const i in filtered_samples) {
        var choiceSelection = document.createElement('input');
        var choiceLabel = document.createElement('label');

        choiceSelection.setAttribute('type', 'radio');
        choiceSelection.setAttribute('name', 'sample_id');
        choiceSelection.setAttribute('id', filtered_samples[i]);

        choiceLabel.innerHTML=filtered_samples[i]+'<br/><br/>';
        choiceLabel.setAttribute('for', 'sample_id');

        sample_list.appendChild(choiceSelection);
        sample_list.appendChild(choiceLabel);
    }

    for (const i in filtered_panels) {
        var choiceSelection = document.createElement('input');
        var choiceLabel = document.createElement('label');

        choiceSelection.setAttribute('type', 'checkbox');
        choiceSelection.setAttribute('name', 'panel_name');
        choiceSelection.setAttribute('id', filtered_panels[i]);

        choiceLabel.innerHTML=filtered_panels[i]+'<br/>';
        choiceLabel.setAttribute('for', 'panel_name');

        panel_list.appendChild(choiceSelection);
        panel_list.appendChild(choiceLabel);
    }

    $("#loader").hide();
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

    var runs_box = document.getElementById("report_runs_box");
    var runs_list = document.getElementById('report_runs_list');

    for (const i in filtered_samples) {
        var choiceSelection = document.createElement('input');
        var choiceLabel = document.createElement('label');

        choiceSelection.setAttribute('type', 'radio');
        choiceSelection.setAttribute('name', 'run_id');
        choiceSelection.setAttribute('id', filtered_samples[i]);

        choiceLabel.innerHTML=filtered_samples[i]+'<br/><br/>';
        choiceLabel.setAttribute('for', 'run_id');

        runs_list.appendChild(choiceSelection);
        runs_list.appendChild(choiceLabel);
    }

    $("#loader").hide();
}