var workflow = '';
var pipeline = '';
var match = '';
var lane = '';
var run_id = '';
var exome = '';
var nextflow_command = '';

function distinct(value, index, self) {
    return self.indexOf(value) === index;
}

function launch_pipeline() {

    $("#command_box").hide();

    console.log(nextflow_command);

    var batch = new AWS.Batch({apiVersion: '2016-08-10'});

    var params = {
        jobDefinition: "nextflow-ufl-germline:10", 
        jobName: "ufl-germline", 
        jobQueue: "hakmonkey-nextflow",
        containerOverrides: {
            'command': [
                'bash',
                '-c',
                nextflow_command+'; aws s3 rm s3://hakmonkey-genetics-lab/Pipieline_Output/_work/'+run_id+'/'
            ]
        }
    };
    batch.submitJob(params, function(err, data) {
        if(err) {
            console.log(err, err.stack);
        } else {
            $("#launch_img").show();
            console.log(data);
        }
    });
}


function check_germline() {

    nextflow_command = 'nextflow run /data/main.nf -work-dir \"s3://hakmonkey-genetics-lab/Pipeline_Output/_work/'+run_id+'/\" --bucket \"s3://hakmonkey-genetics-lab\" --run_id \"'+run_id+'\" --single_lane \"'+lane+'\" --match \"'+match+'\" --exome \"'+exome+'\" --pipeline \"'+pipeline+'\"';

    $("#command_box").show();
    $("#germline_command").show();

    document.getElementById("command_back").setAttribute("onclick", "switchToMatchView()");
    
    if(exome == "NO") {
        var type = "WGS";
    } else if(exome == "YES") {
        var type = "WES";
    }

    document.getElementById('workflow_cell').innerHTML = pipeline;
    document.getElementById('type_cell').innerHTML = type;
    document.getElementById('run_id_cell').innerHTML = run_id;
    document.getElementById('lane_cell').innerHTML = lane;
    document.getElementById('match_cell').innerHTML = match;
}


function check_multiqc() {

    $("#runs_box").hide();
    $("#command_box").show();
    $("#multiqc_command").show();

    document.getElementById("command_back").setAttribute("onclick", "switchToRunIdView()");

    var runs = document.getElementsByName('run_id');
    for(var i in runs) {
        if(runs[i].checked==true){
            run_id = runs[i].id;
        }
    }

    nextflow_command = 'nextflow run /data/main.nf -work-dir \"s3://hakmonkey-genetics-lab/Pipeline_Output/_work/'+run_id+'/\" --bucket \"s3://hakmonkey-genetics-lab\" --run_id \"'+run_id+'\" --pipeline \"'+pipeline+'\"';
    
    document.getElementById('qc_workflow_cell').innerHTML = pipeline;
    document.getElementById('qc_run_id_cell').innerHTML = run_id;
}


function set_match() {

    $("#match_box").hide();

    var match_one = document.getElementById('match_one');
    var match_two = document.getElementById('match_two');
    var other_match = document.getElementById('other_match').value;

    if(match_one.checked==true) {
        match = '_{R1,R2}_001.fastq.gz';

    } else if(match_two.checked==true) {
        match = '_{1,2}.fq.gz';
    } else {
        match = other_match;
    }

    check_germline();
}


function get_match() {

    var single = document.getElementById('single');
    var multi = document.getElementById('multi');

    $("#lane_box").hide();

    if(single.checked==true) {

        lane = 'YES';

        $("#match_box").show();

    } else if(multi.checked==true) {

        lane = 'NO';
        match = '';
        check_germline();
    }
}


function lanes() {

    var runs = document.getElementsByName('run_id');
    for(var i in runs) {
        if(runs[i].checked==true){
            run_id = runs[i].id;
        }
    }

    $("#runs_box").hide();
    $("#lane_box").show();

}

async function get_runs() {

    $("#exome_question_box").hide();

    if(workflow == 'germline_wgs') {
        var prefix = 'Fastqs/';
        var start = 7;
        var stop = 15;
    } else if(workflow == 'germline_wes') {
        var prefix = 'Exome_Fastqs/';
        var start = 13;
        var stop = 21;
    } else if(workflow == 'multiqc') {
        var prefix = 'Pipeline_Output/';
        var start = 16;
        var stop = 24
    }

    s3 = new AWS.S3({apiVersion: '2006-03-01'});

    var bucketParams = {
        Bucket : 'hakmonkey-genetics-lab',
        Prefix: prefix,
        Delimiter: '/'
    };
    
    var samples = [];
    var filtered_samples = [];

    s3.listObjects(bucketParams, function(err, data) {
        if (err) {
            console.log("Error", err);
        } else {
            if(workflow != 'multiqc'){
                for(const i in data['Contents']) {
                    samples.push(data['Contents'][i]['Key'].substring(start,stop));
                }
                samples = samples.filter(distinct);
                samples = samples.filter(item => item);
                for(const i in samples) {
                    filtered_samples.push(samples[i]);
                }
            } else if(workflow == 'multiqc') {
                for(const i in data['CommonPrefixes']) {
                    samples.push(data['CommonPrefixes'][i]['Prefix'].substring(start,stop));
                }
                samples = samples.filter(distinct);
                samples = samples.filter(item => item);
                for(const i in samples) {
                    if(samples[i].startsWith('_') == false) {
                        filtered_samples.push(samples[i]);
                    }
                }
            }
        }
    });

    $("#loader").show();

    await new Promise(r => setTimeout(r, 2000));

    var runs_box = document.getElementById("runs_box");
    var runs_list = document.getElementById('runs_list');

    for (const i in filtered_samples) {
        var choiceSelection = document.createElement('input');
        var choiceLabel = document.createElement('label');

        choiceSelection.setAttribute('type', 'radio');
        choiceSelection.setAttribute('name', 'run_id');
        choiceSelection.setAttribute('id', filtered_samples[i]);

        choiceLabel.appendChild(choiceSelection);
        choiceLabel.innerHTML += '\t'+filtered_samples[i]+'<br/><br/>';
        //choiceLabel.setAttribute('for', 'run_id');

        //runs_list.appendChild(choiceSelection);
        runs_list.appendChild(choiceLabel);
    }

    $("#loader").hide();
    $("#runs_box").show();

    var back_button = document.getElementById("run_box_back");
    var run_button = document.createElement('button');
    var button_label = document.createTextNode('Submit');

    if(workflow != 'multiqc') {
        run_button.setAttribute("onclick", "lanes()");
        back_button.setAttribute("onclick", "switchToTypeView()");
    } else if(workflow == 'multiqc'){
        run_button.setAttribute('onclick', "check_multiqc()");
        back_button.setAttribute("onclick", "switchToPipelineView()");
    }

    run_button.setAttribute('id', 'run_button');
    run_button.appendChild(button_label);
    runs_box.appendChild(run_button);

    
}


function get_type() {
    var WES = document.getElementById("WES");
    var WGS = document.getElementById("WGS");

    $("#exome_question_box").hide;

    if(WGS.checked==true) {

        exome = 'NO';

        workflow = 'germline_wgs';
    } else if(WES.checked==true) {

        exome = 'YES';

        workflow = 'germline_wes';
    }
    get_runs();
}


function get_workflow() {
    var germline = document.getElementById("germline");
    var multiqc = document.getElementById("multiqc");

    $("#workflow_question_box").hide();

    if(germline.checked==true) {

        pipeline = 'GERMLINE';

        $("#exome_question_box").show();

    } else if(multiqc.checked==true) {

        pipeline = 'MULTIQC';
        workflow = 'multiqc';
        get_runs();
    }
}
