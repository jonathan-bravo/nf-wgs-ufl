var run_id = '';

function send_transfer_request() {

    $("#transfer_run_box").hide();

    var batch = new AWS.Batch({apiVersion: '2016-08-10'});

    var params = {
        jobDefinition: "bs-to-aws-ufl-germline:1", 
        jobName: "bs-to-aws_"+run_id, 
        jobQueue: "hakmonkey-bs-to-aws",
        containerOverrides: {
            'command': [
                'bash',
                '-c',
                './bs-to-aws.sh -b hakmonkey-genetics-lab -r '+run_id+';'
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

    //$("#menu").show();
    //location.reload();
}

function request_run_transfer() {

    run_id = document.getElementById('transfer_run_id').value;

    send_transfer_request();
}
