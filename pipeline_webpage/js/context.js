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
}

function switchToDocumentationView(){
    window.location.replace("toc.html")
}