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
import { storageGet, storageSet, KEYS } from "../util/localStorage";

@connect((state) => ({
  tosURL: state.config?.parameters?.about_legal_tos,
  privacyURL: state.config?.parameters?.about_legal_privacy,
}))
class TermsPrompt extends React.PureComponent {
  constructor(props) {
    super(props);
    const { tosURL, privacyURL } = this.props;
    const cookieDecision = storageGet(KEYS.COOKIE_DECISION, null);
    const hasDecided = cookieDecision !== null;
    this.state = {
      hasDecided,
      isEnabled: !!tosURL || !!privacyURL,
      isOpen: !hasDecided,
    };
  }

  componentDidMount() {
    const { hasDecided, isEnabled } = this.state;
    if (isEnabled && !hasDecided) {
      this.setState({ isOpen: true });
    }
  }

  handleOK = () => {
    this.setState({ isOpen: false });
    storageSet(KEYS.COOKIE_DECISION, "yes");
    if (window.cookieDecisionCallback instanceof Function) {
      try {
        window.cookieDecisionCallback();
      } catch (e) {
        // continue
      }
    }
  };

  handleNo = () => {
    this.setState({ isOpen: false });
    storageSet(KEYS.COOKIE_DECISION, "no");
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
          rel="noopener"
        >
          terms of service (updates coming)
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
          rel="noopener"
        >
          privacy policy (updates coming)
        </a>
        .&nbsp;
      </span>
    );
  }

  render() {
    const { isOpen, isEnabled } = this.state;
    if (!isEnabled || !isOpen) return null;
    return (
      <Drawer
        onclose={this.drawerClose}
        isOpen={isOpen}
        size="120px"
        position={Position.BOTTOM}
        canOutsideClickClose={false}
        hasBackdrop /* if the user can't use the app or click outside to dismiss, that should be visually represented with a backdrop */
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
            <Button
              intent="primary"
              onClick={this.handleOK}
              data-testid="tos-cookies-accept"
            >
              I&apos;m OK with cookies!
            </Button>{" "}
            <Button onClick={this.handleNo} data-testid="tos-cookies-reject">
              No thanks
            </Button>
          </div>
        </div>
      </Drawer>
    );
  }
}

export default TermsPrompt;
