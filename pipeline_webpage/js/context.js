function switchToLogInView(){
    $("#userNameInput").val('');
    $("#passwordInput").val('');
    $("#userNameInput").show();
    $("#passwordInput").show();
    $("#logInButton").show();
    $("#logOutButton").hide();
    $("#pipeline").hide();
    $("#documentation").hide();
    $("#workflow_question_box").hide();
    $("#exome_question_box").hide();
    $("#runs_box").hide();
    $("#lane_box").hide();
    $("#match_box").hide();
    $("#command_box").hide();
    $("#launch_img").hide();
    $("#report").hide();
    $("#report_runs_box").hide();
}

function switchToLoggedInView(){
    $("#userNameInput").hide();
    $("#passwordInput").hide();
    $("#logInButton").hide();
    $("#logOutButton").show();
    $("#pipeline").show();
    $("#documentation").show();
    $("#workflow_question_box").hide();
    $("#exome_question_box").hide();
    $("#runs_box").hide();
    $("#lane_box").hide();
    $("#match_box").hide();
    $("#command_box").hide();
    $("#launch_img").hide();
    $("#report").show();
    $("#report_runs_box").hide();
}

function switchToPipelineView(){
    $("#userNameInput").hide();
    $("#passwordInput").hide();
    $("#logInButton").hide();
    $("#logOutButton").show();
    $("#pipeline").hide();
    $("#documentation").hide();
    $("#workflow_question_box").show();
    $("#exome_question_box").hide();
    $("#runs_box").hide();
    $("#lane_box").hide();
    $("#match_box").hide();
    $("#command_box").hide();
    $("#launch_img").hide();
    $("#report").hide();
    $("#report_runs_box").hide();
}

function switchToReportView(){
    $("#userNameInput").hide();
    $("#passwordInput").hide();
    $("#logInButton").hide();
    $("#logOutButton").show();
    $("#pipeline").hide();
    $("#documentation").hide();
    $("#workflow_question_box").hide();
    $("#exome_question_box").hide();
    $("#runs_box").hide();
    $("#lane_box").hide();
    $("#match_box").hide();
    $("#command_box").hide();
    $("#launch_img").hide();
    $("#report").hide();
    $("#report_runs_box").show();
    get_report_runs();
}

function switchToDocumentationView(){
    window.location.replace("toc.html")
}