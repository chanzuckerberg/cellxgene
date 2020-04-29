import React from "react";
import { connect } from "react-redux";
import * as globals from "../../globals";
import { termsOfServiceToast } from "../framework/toasters";

const TosDismissedKey = "cxg.tosDismissed";

function storageGet(key, defaultValue = null) {
  try {
    const val = window.localStorage.getItem(key);
    if (val === null) return defaultValue;
    return val;
  } catch (e) {
    return defaultValue;
  }
}

function storageSet(key, value) {
  try {
    window.localStorage.setItem(key, value);
  } catch {
    
  }
}

@connect(state => ({
  tosURL: state.config?.parameters?.about_legal_tos
}))
class TermsPrompt extends React.PureComponent {
  constructor(props) {
    super(props);
    this.state = {
      hasDismissed: storageGet(TosDismissedKey, false)
    };
  }

  componentDidMount() {
    const { hasDismissed } = this.state;
    const { tosURL } = this.props;
    if (!hasDismissed && tosURL) {
      this.popTermsToast();
    }
  }

  onTermsToastDismissed = () => {
    this.setState({ hasDismissed: "yes" });
    storageSet(TosDismissedKey, "yes");
  };

  popTermsToast() {
    const { tosURL } = this.props;
    termsOfServiceToast(
      <span>
        By using our site, you are agreeing to our{" "}
        <a
          style={{
            fontWeight: 700,
            color: "white",
            textDecoration: "underline"
          }}
          href={tosURL}
          target="_blank"
        >
          Terms of Service
        </a>
      </span>,
      this.onTermsToastDismissed
    );
  }

  render() {
    return null;
  }
}

export default TermsPrompt;
