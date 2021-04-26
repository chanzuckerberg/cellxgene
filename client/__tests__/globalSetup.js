const {
  SecretsManagerClient,
  GetSecretValueCommand,
} = require("@aws-sdk/client-secrets-manager");

const { setup } = require("jest-environment-puppeteer");

const client = new SecretsManagerClient({ region: "us-west-2" });

const secretValueRequest = {
  SecretId: "corpora/backend/dev/auth0-secret",
};

const command = new GetSecretValueCommand(secretValueRequest);

module.exports = async () => {
  await setup();
  try {
    const secret = JSON.parse((await client.send(command)).SecretString);
    process.env.TEST_ACCOUNT_PASS = secret.test_account_password;
  } catch (error) {
    console.error(error);
  }
  console.log(process.env.TEST_ACCOUNT_PASS);
};
