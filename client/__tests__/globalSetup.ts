const {
  SecretsManagerClient,
  GetSecretValueCommand,
  // eslint-disable-next-line @typescript-eslint/no-var-requires --- FIXME: disabled temporarily on migrate to TS.
} = require("@aws-sdk/client-secrets-manager");

// eslint-disable-next-line @typescript-eslint/no-var-requires --- FIXME: disabled temporarily on migrate to TS.
const { setup } = require("jest-environment-puppeteer");

const client = new SecretsManagerClient({ region: "us-west-2" });

const deploymentStage = process.env.DEPLOYMENT_STAGE || "test";

const secretValueRequest = {
  SecretId: `corpora/backend/${deploymentStage}/auth0-secret`,
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
};
