import React from "react";
import { connect } from "react-redux";
import {
  Drawer,
  Button,
  Classes,
  Position,
  Colors,
  Icon,
} from "@blueprintjs/core";
import * as globals from "../../globals";

const CookieDecision = "cxg.cookieDecision";

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
    return;
  }
}

@connect((state) => ({
  tosURL: state.config?.parameters?.about_legal_tos,
  privacyURL: state.config?.parameters?.about_legal_privacy,
}))
class TermsPrompt extends React.PureComponent {
  constructor(props) {
    super(props);
    const cookieDecision = storageGet(CookieDecision, null);
    const hasDecided = cookieDecision !== null;
    this.state = {
      hasDecided,
      cookieDecision,
      isOpen: !hasDecided,
    };
  }

  componentDidMount() {
    const { hasDecided } = this.state;
    const { tosURL, privacyURL } = this.props;
    if (!hasDecided && (tosURL || privacyURL)) {
      this.openConsentDrawer();
    }
  }

  openConsentDrawer() {
    this.setState({ isOpen: true });
  }

  closeConsentDrawer() {
    this.setState({ isOpen: false });
  }

  handleOK = () => {
    this.setState({ isOpen: false });
    storageSet(CookieDecision, "yes");
  };

  handleNo = () => {
    this.setState({ isOpen: false });
    storageSet(CookieDecision, "no");
  };

  renderTos() {
    const { tosURL } = this.props;
    if (!tosURL) return null;
    return (
      <span>
        <Icon icon="info-sign" intent="primary" /> By using this site, you are
        agreeing to our{" "}
        <a
          style={{
            fontWeight: 700,
            textDecoration: "underline",
          }}
          href={tosURL}
          target="_blank"
        >
          terms of service
        </a>
        .{" "}
      </span>
    );
  }

  renderPrivacy() {
    const { privacyURL } = this.props;
    if (!privacyURL) return null;
    return (
      <span>
        To learn more, read our{" "}
        <a
          style={{
            fontWeight: 700,
            textDecoration: "underline",
          }}
          href={privacyURL}
          target="_blank"
        >
          privacy policy
        </a>
        .&nbsp;
      </span>
    );
  }

  render() {
    const { tosURL, privacyURL } = this.props;
    const { isOpen } = this.state;
    return (
      <Drawer
        onclose={this.drawerClose}
        isOpen={isOpen}
        size={"120px"}
        position={Position.BOTTOM}
        canOutsideClickClose={false}
        hasBackdrop={false}
        enforceFocus={false}
        autoFocus={false}
        portal={false}
      >
        <div
          className={Classes.DRAWER_BODY}
          style={{ backgroundColor: Colors.LIGHT_GRAY1 }}
        >
          <div className={Classes.DIALOG_BODY}>
            <div>
              {this.renderTos()}
              <span>
                We use cookies to help us improve the site and to inform our
                future efforts, and we also use necessary cookies to make our
                site work.&nbsp;
              </span>
              {this.renderPrivacy()}
            </div>
          </div>
          <div className={Classes.DIALOG_FOOTER} style={{ textAlign: "left" }}>
            <Button intent="primary" onClick={this.handleOK}>
              I'm OK with cookies!
            </Button>{" "}
            <Button onClick={this.handleNo}>No thanks</Button>
          </div>
        </div>
      </Drawer>
    );
  }
}

export default TermsPrompt;
