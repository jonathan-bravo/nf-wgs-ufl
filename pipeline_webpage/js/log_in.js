//=============== AWS IDs ===============
var userPoolId = 'us-east-1_yu5nvjI1j';
var clientId = '52of22i4vf76flprhik1h5rn7';
var region = 'us-east-1';
var identityPoolId = 'us-east-1:9dd82334-7556-41ab-b62e-dd6fe51a0134';
//=============== AWS IDs ===============

var cognitoUser;
var idToken;
var userPool;

var poolData = { 
    UserPoolId : userPoolId,
    ClientId : clientId
};

getCurrentLoggedInSession();

function logOut(){
    if (cognitoUser != null) {

        $("#loader").show();
        cognitoUser.signOut();
        switchToLogInView();
        console.log('Logged out!');
        $("#loader").hide();
        location.reload();
    }
}

/*
Starting point for user login flow with input validation
*/
function logIn(){

    if(!$('#userNameInput').val() || !$('#passwordInput').val()){
        console.log('Please enter Username and Password!');
    }else{
        var authenticationData = {
            Username : $('#userNameInput').val(),
            Password : $("#passwordInput").val(),
        };
        var authenticationDetails = new AmazonCognitoIdentity.AuthenticationDetails(authenticationData);

        var userData = {
            Username : $('#userNameInput').val(),
            Pool : userPool
        };
        cognitoUser = new AmazonCognitoIdentity.CognitoUser(userData);

        $("#loader").show();
        cognitoUser.authenticateUser(authenticationDetails, {
            onSuccess: function (result) {
                console.log('Logged in!');
                switchToLoggedInView();

                idToken = result.getIdToken().getJwtToken();
                getCognitoIdentityCredentials();
            },

            onFailure: function(err) {
                console.log(err.message);
                $("#loader").hide();
            },

        });
    }
}

/*
This method will get temporary credentials for AWS using the IdentityPoolId and the Id Token recieved from AWS Cognito authentication provider.
*/
function getCognitoIdentityCredentials(){
    AWS.config.region = region;

    var loginMap = {};
    loginMap['cognito-idp.' + region + '.amazonaws.com/' + userPoolId] = idToken;

    AWS.config.credentials = new AWS.CognitoIdentityCredentials({
        IdentityPoolId: identityPoolId,
        Logins: loginMap
    });

    AWS.config.credentials.clearCachedId();

    AWS.config.credentials.get(function(err) {
        if (err){
            console.log(err.message);
        }
        else {
            console.log('AWS Access Key: '+ AWS.config.credentials.accessKeyId);
            console.log('AWS Secret Key: '+ AWS.config.credentials.secretAccessKey);
            console.log('AWS Session Token: '+ AWS.config.credentials.sessionToken);
        }

        $("#loader").hide();
    });
}

/*
If user has logged in before, get the previous session so user doesn't need to log in again.
*/
function getCurrentLoggedInSession(){

    $("#loader").show();
    userPool = new AmazonCognitoIdentity.CognitoUserPool(poolData);
    cognitoUser = userPool.getCurrentUser();

    if(cognitoUser != null){
        cognitoUser.getSession(function(err, session) {
            if (err) {
                console.log(err.message);
            }else{
                console.log('Session found! Logged in.');
                switchToLoggedInView();
                idToken = session.getIdToken().getJwtToken();
                getCognitoIdentityCredentials();
            }
            $("#loader").hide();
        });
    }else{
        console.log('Session expired. Please log in again.');
        $("#loader").hide();
    }

}