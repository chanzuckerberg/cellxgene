import React from "react";
import { connect } from "react-redux";
import * as globals from "../../globals";
import {
  termsOfServiceToast
} from "../framework/toasters";

@connect(state => ({
  tosURL: state.config?.parameters?.about_legal_tos,
}))
class TermsPrompt extends React.PureComponent {
  constructor(props) {
    super(props);
    this.state = {
    };
  }

  componentDidMount () {
   if (/* user has not agreed... */ && tosURL) { 
     this.popTermsToast() 
    }
  }

  onTermsToastDismissed () {
    /* Bruce to fire Action */
  }

  popTermsToast () {
    termsOfServiceToast(
      (<span>
        By using our site, you are agreeing to our{" "} 
        <a style={{
          fontWeight: 700, 
          color: "white", 
          textDecoration: "underline"
        }} 
        href={tosURL}>
          Terms of Service
        </a>
      </span>), this.onTermsToastDismissed
    )
  }

  render() {
    return null

  }
}

export default TermsPrompt;
